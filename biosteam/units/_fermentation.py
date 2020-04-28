# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:45:47 2018

@author: yoelr
"""
import numpy as np
from ._hx import HXutility
from .. import Unit
from scipy.integrate import odeint
from .decorators import cost
from .design_tools import size_batch
from thermosteam.reaction import Reaction

__all__ = ('Fermentation',)

@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3, S=3785, n=0.6, BM=1.8, N='N')
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='N')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000,
      S=3785, n=0.5, BM=1.5, N='N')
class Fermentation(Unit):
    """
    Create a Fermentation object which models large-scale batch fermentation
    for the production of 1st generation ethanol using yeast
    [1]_ [2]_ [3]_ [4]_. A compound with CAS 'Yeast' must be present.
    Only sucrose and glucose are taken into account for conversion.
    Conversion is based on reaction time, `tau`. Cleaning and unloading time,
    `tau_0`, fraction of working volume, `V_wf`, and number of reactors,
    `N_reactors`, are attributes that can be changed. Cost of a reactor
    is based on the NREL batch fermentation tank cost assuming volumetric
    scaling with a 6/10th exponent [5]_. 
    
    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed into the fermentor.
    outs : stream sequence
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time.
    efficiency=0.9 : float, optional
        User enforced efficiency.
    iskinetic=False: bool, optional
        If True, `Fermenation.kinetic_model` will be used.
    N : int
        Number of batch reactors
    T=305.15 : float
        Temperature of reactor [K].
    
    Examples
    --------
    Simulate a Fermentation object which models batch fermentation for the
    production of 1st generation ethanol using yeast.
    
    >>> from biorefineries.lipidcane.chemicals import ethanol_chemicals 
    >>> from biosteam.units import Fermentation
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(ethanol_chemicals)
    >>> feed = Stream('feed',
    ...               Water=1.20e+05,
    ...               Glucose=1.89e+03,
    ...               Sucrose=2.14e+04,
    ...               DryYeast=1.03e+04,
    ...               units='kg/hr',
    ...               T=32+273.15)
    >>> F1 = Fermentation('F1',
    ...                   ins=feed, outs=('CO2', 'product'),
    ...                   tau=8, efficiency=0.90, N=8)
    >>> F1.simulate()
    >>> F1.show()
    Fermentation: F1
    ins...
    [0] feed
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water     6.66e+03
                        Glucose   10.5
                        Sucrose   62.5
                        DryYeast  1.03e+04
    [1] missing stream
    outs...
    [0] CO2
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    2.5
                        CO2      244
                        Ethanol  0.582
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water     6.6e+03
                        Ethanol   243
                        Glucose   13.6
                        DryYeast  1.03e+04
    >>> F1.results()
    Fermentation                                       Units        F1
    Power               Rate                              kW      11.6
                        Cost                          USD/hr     0.908
    Chilled water       Duty                           kJ/hr -2.91e+06
                        Flow                         kmol/hr  1.95e+03
                        Cost                          USD/hr      14.5
    Design              Reactor volume                    m3       246
                        Batch time                        hr      12.6
                        Loading time                      hr      1.57
                        Cleaning and unloading time       hr         3
                        Working volume fraction                    0.9
                        Number of reactors                           8
    Purchase cost       Coolers                          USD  1.62e+05
                        Reactors                         USD  1.87e+06
                        Agitators                        USD  1.16e+05
                        Cleaning in place                USD   7.1e+05
    Total purchase cost                                  USD  2.86e+06
    Utility cost                                      USD/hr      15.4
    
    References
    ----------
    .. [1] Oliveira, Samuel C., et al. "Discrimination between ethanol inhibition models in a continuous alcoholic fermentation process using flocculating yeast." Applied biochemistry and biotechnology 74.3 (1998): 161-172.
    
    .. [2] Oliveira, Samuel C., et al. "Continuous ethanol fermentation in a tower reactor with flocculating yeast recycle: scale-up effects on process performance, kinetic parameters and model predictions." Bioprocess Engineering 20.6 (1999): 525-530.
    
    .. [3] Oliveira, Samuel C., et al. "Mathematical modeling of a continuous alcoholic fermentation process in a two-stage tower reactor cascade with flocculating yeast recycle." Bioprocess and biosystems engineering 38.3 (2015): 469-479.
    
    .. [4] Oliveira, Samuel C., et al. "Kinetic Modeling of 1‐G Ethanol Fermentations." Fermentation Processes. InTech, 2017.
    
    .. [5] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden National. Renewable Energy Laboratory Golden, Colorado. P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical Report NREL/TP-5100-47764
    
    """
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Batch time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr'}
    _N_ins = _N_outs = 2
    _N_heat_utilities = 1
    _has_power_utility = True
    line = 'Fermentation'    
    
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr)
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume
    V_wf = 0.9
    
    #: tuple[float] Kinetic parameters for the kinetic model. Default constants are fitted for Oliveria's model (mu_m1, mu_m2, Ks1, Ks2, Pm1, Pm2, Xm, Y_PS, a)
    kinetic_constants = (0.31,  # mu_m1
                         1.01,  # mu_m2
                         1.88,  # Ks1
                         2.81,  # Ks2
                         82.8,  # Pm1
                         108.2, # Pm2
                         113.4, # Xm
                         0.45,  # Y_PS
                         0.18)  # a
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''),
                ('Number of reactors', self.N, ''))
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 tau,  N, efficiency=0.9, iskinetic=False, T=305.15):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.hydrolysis = Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00)
        self.fermentation = Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', efficiency)
        self.iskinetic = iskinetic
        self.efficiency = efficiency
        self.tau = tau
        self.N = N
        self.T = T
        self.cooler = HXutility(None)
        
    def _setup(self):
        vent, effluent = self.outs
        vent.phase = 'g'
        self.cooler._outs[0].T = effluent.T = vent.T = self.T

    def _calc_efficiency(self, feed, tau):
        # Get initial concentrations
        y, e, s, w = feed.indices(['Yeast',
                                   '64-17-5',
                                   '492-61-5',
                                   '7732-18-5'])
        mass = feed.mass
        F_vol = feed.F_vol
        concentration_in = mass/F_vol
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
    def kinetic_model(z, t, *kinetic_constants):
        """
        Return change of yeast, ethanol, and substrate concentration in kg/m3.
        
        Parameters
        ----------
        z : Iterable with (X, E, S) [-]:
            * X: Yeast concentration (kg/m3)
            * P: Ethanol concentration (kg/m3)
            * S: Substrate concentration (kg/m3)
        
        t : float
            Time point
        
        *kinetic_constants
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
        effluent.mix_from(self.ins)
        effluent_mol = effluent.mol
        self.hydrolysis(effluent_mol)
        if self.iskinetic:
            self.fermentation.X = self._calc_efficiency(effluent, self._tau)
        self.fermentation(effluent_mol)
        vent.receive_vent(effluent)
    
    @property
    def N_at_minimum_capital_cost(self):
        cost_old = np.inf
        self.autoselect_N = False
        self._N, N = 2, self._N
        cost_new = self.purchase_cost
        self._summary()
        while cost_new < cost_old:
            self._N += 1
            self._summary()
            cost_old = cost_new
            cost_new = self.purchase_cost
        self._N, N = N, self._N
        self.autoselect_N = True
        return N - 1
        
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        tau = self._tau
        tau_0 = self.tau_0
        Design = self.design_results
        if self.autoselect_N:
            self._N = self.N_at_minimum_capital_cost
        N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, self.V_wf))
        hx_effluent = effluent.copy()
        hx_effluent.mol[:] /= N
        cooler = self.cooler
        cooler.simulate_as_auxiliary_exchanger(self.Hnet/N, hx_effluent)
        hu_fermentation, = self.heat_utilities
        hu_cooler, = cooler.heat_utilities
        hu_fermentation.copy_like(hu_cooler)
        hu_fermentation.scale(N)
        self.purchase_costs['Coolers'] = self.cooler.purchase_costs['Heat exchanger'] * N
        
    