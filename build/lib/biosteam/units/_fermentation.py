# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:45:47 2018

@author: yoelr
"""
import numpy as np
from ._hx import HXutility
from .. import Stream, Unit
from scipy.integrate import odeint
from .decorators import cost
from ._tank import MixTank
from .designtools import size_batch
from ..reaction import Reaction

@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3, S=3785, n=0.6, BM=1.8)
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500,
      S=3785, n=0.5, kW=22.371, BM=1.5)
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000,
      S=3785, n=0.5, BM=1.5)
class Fermentation(Unit):
    """Create a Fermentation object which models large-scale batch fermentation for the production of 1st generation ethanol using yeast [1, 2, 3, 4]. Only sucrose and glucose are taken into account. Conversion is based on reaction time, `tau`. Cleaning and unloading time, `tau_0`, fraction of working volume, `V_wf`, and number of reactors, `N_reactors`, are attributes that can be changed. Cost of a reactor is based on the NREL batch fermentation tank cost assuming volumetric scaling with a 6/10th exponent [3].
    
    **Parameters**
    
        **tau:** Reaction time.
        
        **efficiency:** User enforced efficiency.
        
        **N:** Number of batch reactors
    
    **ins**
    
        [:] Inffluent streams
        
    **outs**
    
        [0] Effluent
        
        [1] CO2
    
    **References**
    
        [1] Oliveira, Samuel C., et al. "Discrimination between ethanol inhibition models in a continuous alcoholic fermentation process using flocculating yeast." Applied biochemistry and biotechnology 74.3 (1998): 161-172.
        
        [2] Oliveira, Samuel C., et al. "Continuous ethanol fermentation in a tower reactor with flocculating yeast recycle: scale-up effects on process performance, kinetic parameters and model predictions." Bioprocess Engineering 20.6 (1999): 525-530.
        
        [3] Oliveira, Samuel C., et al. "Mathematical modeling of a continuous alcoholic fermentation process in a two-stage tower reactor cascade with flocculating yeast recycle." Bioprocess and biosystems engineering 38.3 (2015): 469-479.
        
        [4] Oliveira, Samuel C., et al. "Kinetic Modeling of 1‐G Ethanol Fermentations." Fermentation Processes. InTech, 2017.
        
        [5] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden National. Renewable Energy Laboratory Golden, Colorado. P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical Report NREL/TP-5100-47764
        
    .. note::
        
        A compound with CAS 'Yeast' must be present in species.
        
    **Examples**
    
        :doc:`Fermentation Example`
        
    """
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr'}
    _N_ins = _N_outs = 2
    _N_heat_utilities = 0
    _has_power_utility = True
    line = 'Fermentation'    
    
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: Cleaning and unloading time (hr)
    tau_0 = 3
    
    #: Fraction of filled tank to total tank volume
    working_volume_fraction = MixTank.working_volume_fraction
    _V_wf = 0.9
    
    #: tuple of kinetic parameters for the kinetic model. Default constants are fitted for Oliveria's model (mu_m1, mu_m2, Ks1, Ks2, Pm1, Pm2, Xm, Y_PS, a)
    kinetic_constants = (0.31,  # mu_m1
                         1.01,  # mu_m2
                         1.88,  # Ks1
                         2.81,  # Ks2
                         82.8,  # Pm1
                         108.2, # Pm2
                         113.4, # Xm
                         0.45,  # Y_PS
                         0.18)  # a
    
    def __init__(self, ID='', ins=None, outs=(), *, 
                 tau,  N, efficiency=0.9, iskinetic=False):
        Unit.__init__(self, ID, ins, outs)
        self.hydrolysis = Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00)
        self.fermentation = Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', efficiency)
        self.iskinetic = iskinetic
        self.efficiency = efficiency
        self.tau = tau
        self.N = N
        self._cooler = hx = HXutility(None)
        self._heat_utilities = hx._heat_utilities
        hx._ins = hx._outs
        vent, effluent = self.outs
        hx._outs[0].T = effluent.T = vent.T = 305.15
        vent.phase = 'g'

    def _calc_efficiency(self, feed, tau):
        # Get initial concentrations
        y, e, s, w = feed.indices(['Yeast',
                                   '64-17-5',
                                   '492-61-5',
                                   '7732-18-5'])
        mass = feed.mass
        volnet = feed.volnet
        concentration_in = mass/volnet
        X0, P0, S0 = (concentration_in[i] for i in (y, e, s))
        
        # Integrate to get final concentration
        t = np.linspace(0, tau, 1000)
        C_t = odeint(self.kinetic_model, (X0, P0, S0), t,
                     args=self.kinetic_constants)
        # Cache data
        self._X = C_t[:, 0]
        self._P = C_t[:, 1]
        self._S = S = C_t[:, 2]
        
        # Calculate efficiency
        Sf = S[-1]
        Sf = Sf if Sf > 0 else 0
        Y_PS = self.kinetic_constants[-2]
        eff = (S0 - Sf)/S0 * Y_PS/0.511
        return eff
        
    @staticmethod
    def kinetic_model(z, t, *kinetic_constants) -> '(dXdt, dPdt, dSdt)':
        """Return change of yeast, ethanol, and substrate concentration in kg/m3.
        
        **Parameters**
        
            **z:** Iterable with (X, E, S) [-]:
                * X: Yeast concentration (kg/m3)
                * P: Ethanol concentration (kg/m3)
                * S: Substrate concentration (kg/m3)
            
            **t:** Time point
            
            ***kinetic_constants:**
                * mu_m1: Maximum specific growth rate (1/hr)
                * mu_m2: Maximum specific ethanol production rate (g-product/g-cell-hr)
                * Ks1: Sugar saturation constant for growth (g/L)
                * Ks2: Sugar saturation constant for product (g/L)
                * Pm1: Maximum product concentration at zero growth [mu_m1=0] (g/L)
                * Pm2: Maximum product concentration [mu_m2=0] (g/L)
                * Xm: Maximum cell concentration [mu_m1=0] (g/L)
                * Y_PS: Ethanol yield based on sugar consumed
                * a: Toxic power
                
        """
        mu_m1, mu_m2, Ks1, Ks2, Pm1, Pm2, Xm, Y_PS, a = kinetic_constants
        
        # Current yeast, ethanol, and glucose concentration (kg/m3)
        X, P, S = z
        
        # Compute coefficients
        mu_X = mu_m1 * (S/(Ks1 + S)) * (1 - P/Pm1)**a*((1-X/Xm))
        mu_P = mu_m2 * (S/(Ks2 + S)) * (1 - P/Pm2)
        mu_S = mu_P/0.45
        
        # Compute derivatives
        dXdt = mu_X * X
        dPdt = (mu_P * X)
        dSdt =  - mu_S * X
        return (dXdt, dPdt, dSdt)
       
    @property
    def N(self):
        """[int] Number of reactors"""
        return self._N
    @N.setter
    def N(self, N):
        if N <= 1:
            raise ValueError(f"number of reactors must be greater than 1, value {N} is infeasible")
        self._N = N

    @property
    def efficiency(self):
        return self.fermentation.X
    @efficiency.setter
    def efficiency(self, efficiency):
        self.fermentation.X = efficiency

    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau

    def _run(self):
        vent, effluent = self.outs
        Stream.sum(effluent, self.ins)
        effluent_mol = effluent.mol
        self.hydrolysis(effluent_mol)
        if self.iskinetic:
            self.fermentation.X = self._calc_efficiency(effluent, self._tau)
        self.fermentation(effluent_mol)
        vent.copyflow(effluent, ('CO2',), remove=True)
        vent.recieve_vent(effluent)
    
    @property
    def N_at_minimum_capital_cost(self):
        cost_old = np.inf
        self._N, N = 2, self._N
        cost_new = self.purchase_cost
        self._summary()
        while cost_new < cost_old:
            self._N += 1
            self._summary()
            cost_old = cost_new
            cost_new = self.purchase_cost
        self._N, N = N, self._N
        return N - 1
        
    def _design(self):
        v_0 = self.outs[1].volnet
        tau = self._tau
        tau_0 = self.tau_0
        Design = self._Design
        if self.autoselect_N:
            self.autoselect_N = False
            self._N = self.N_at_minimum_capital_cost
            self.autoselect_N = True
        Design['Number of reactors'] = N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, self._V_wf))
        hx = self._cooler
        hx.outs[0]._mol[:] = self.outs[0].mol/N 
        hu = hx._heat_utilities[0]
        hu(self._Hnet/N, self.outs[0].T)
        hx._design(hu.duty)
        hx._cost()
        hu.duty *= N
        hu.cost *= N
        hu.flow *= N
        
    def _end_decorated_cost_(self):
        N = self._Design['Number of reactors']
        Cost = self._Cost
        Cost['Coolers'] = self._cooler._Cost['Heat exchanger']
        for i in Cost: Cost[i] *= N
    