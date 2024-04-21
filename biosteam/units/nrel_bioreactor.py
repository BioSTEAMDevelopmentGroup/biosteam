# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.nrel_bioreactor.NRELBatchBioreactor
.. autoclass:: biosteam.units.nrel_bioreactor.NRELFermentation

References
----------
.. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
    National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
    J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
    Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
    for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
    Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
    Report NREL/TP-5100-47764

.. [2] Oliveira, Samuel C., et al. "Discrimination between ethanol 
    inhibition models in a continuous alcoholic fermentation process using
    flocculating yeast." Applied biochemistry and biotechnology 74.3 (1998): 161-172.

.. [3] Oliveira, Samuel C., et al. "Continuous ethanol fermentation in a
    tower reactor with flocculating yeast recycle: scale-up effects on process
    performance, kinetic parameters and model predictions." Bioprocess
    Engineering 20.6 (1999): 525-530.

.. [4] Oliveira, Samuel C., et al. "Mathematical modeling of a continuous
    alcoholic fermentation process in a two-stage tower reactor cascade with
    flocculating yeast recycle." Bioprocess and biosystems engineering 38.3
    (2015): 469-479.

.. [5] Oliveira, Samuel C., et al. "Kinetic Modeling of 1‐G Ethanol
    Fermentations." Fermentation Processes. InTech, 2017.
    
    
"""
import numpy as np
from .. import Unit
from .design_tools import size_batch
from .decorators import cost
from math import ceil
from scipy.integrate import odeint
from thermosteam.reaction import Reaction, ParallelReaction

__all__ = (
    'NRELBatchBioreactor', 'NRELFermentation',
    'BatchBioreactor', 'Fermentation', # For backwards compatibility
) 

# %% NREL

@cost('Recirculation flow rate', 'Recirculation pumps', kW=30, S=77.22216,
      cost=47200, n=0.8, BM=2.3, CE=522, N='Number of reactors')
@cost('Reactor volume', 'Cleaning in place', CE=521.9,
      cost=421e3, S=3785, n=0.6, BM=1.8)
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500,
      S=3785, n=0.5, kW=22.371, BM=1.5, N='Number of reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000,
      S=3785, n=0.5, BM=1.5, N='Number of reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=20920000.0, n=0.7, BM=2.2, N='Number of reactors',
      magnitude=True) # Based on a similar heat exchanger
class NRELBatchBioreactor(Unit, isabstract=True):
    """
    Abstract Bioreactor class. Conversion is based on reaction time, `tau`.
    Cleaning and unloading time,`tau_0`, fraction of working volume, `V_wf`,
    and number of reactors, `N_reactors`, are attributes that can be changed.
    The cost of a reactor is based on the NREL batch fermentation tank cost 
    assuming volumetric scaling with a 6/10th exponent [1]_. 
    
    **Abstract methods**
    
    _run():
        Must set outlet stream flows.
        
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed into the reactor.
    outs : 
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time [hr].
    N : int, optional
        Number of batch reactors
    V : float, optional
        Target volume of reactors [m^3].
    T=305.15 : float
        Operating temperature of reactor [K].
    P=101325 : float
        Operating pressure of reactor [Pa].
    Nmin=2 : int
        Minimum number of fermentors.
    Nmax=36: int
        Maximum number of fermentors.  
    
    Notes
    -----
    Either N or V must be given.
    
    """
    _units = {'Reactor volume': 'm3',
              'Cycle time': 'hr',
              'Batch time': 'hr',
              'Loading time': 'hr',
              'Total dead time': 'hr',
              'Reactor duty': 'kJ/hr',
              'Recirculation flow rate': 'm3/hr'}
    _N_ins = _N_outs = 2
    
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, ''))
    
    def _init(self, tau=None, N=None, V=None, T=305.15, P=101325,
              Nmin=2, Nmax=36):
        self._N = N; self._V = V
        
        #: [float] Reaction time [hr].
        self.tau = tau
        
        #: [int] Number of batch reactors
        if N: self.N = N
        
        #: [float] Target volume of a fermentor
        if V: self.V = V
        
        #: [float] Operating temperature of reactor [K].
        self.T = T
        
        #: [float] Operating pressure of reactor [Pa].
        self.P = P
        
        #: [int] Minimum number of fermentors
        self.Nmin = Nmin
        
        #: [int] Maximum number of fermentors
        self.Nmax = Nmax
        
    def _setup(self):
        super()._setup()
        vent, effluent = self.outs
        vent.phase = 'g'
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P
        
    @property
    def vent(self):
        return self.outs[0]
    @property
    def effluent(self):
        return self.outs[1]
        
    @property
    def N(self):
        """[int] Number of reactor."""
        return self._N
    @N.setter
    def N(self, N):
        if N is None: 
            self._N = N
            return
        if N <= 1:
            raise ValueError(f"number of reactors must be greater than 1, value {N} is infeasible")
        assert not self._V, 'cannot specify both reactor volume and number of reactors'
        self._N = ceil(N)

    @property
    def V(self):
        """[float] Reactor volume."""
        return self._V
    @V.setter
    def V(self, V):
        if V <= 1:
            raise ValueError(f"reactor volume must be greater than 1, value {V} is infeasible")
        assert not self._N, 'cannot specify both reactor volume and number of reactors'
        self._V = V

    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, tau):
        self._tau = tau
    
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
        effluent = self.effluent
        v_0 = effluent.F_vol
        tau = self._tau
        tau_0 = self.tau_0
        V_wf = self.V_wf
        Design = self.design_results
        if self.autoselect_N:
            N = self.N_at_minimum_capital_cost
        elif self.V:
            N = v_0 / self.V / V_wf * (tau + tau_0) + 1
            if N < 2:
                N = 2
            else:
                N = ceil(N)
        else:
            N = self._N
        Design.update(size_batch(v_0, tau, tau_0, N, V_wf))
        Design['Number of reactors'] = N
        Design['Recirculation flow rate'] = v_0 / N
        duty = self.Hnet
        Design['Reactor duty'] = duty
        self.add_heat_utility(duty, self.T)


class NRELFermentation(NRELBatchBioreactor):
    """
    Create a Fermentation object which models large-scale batch fermentation
    for the production of 1st generation ethanol using yeast
    [2]_ [3]_ [4]_ [5]_. A compound with CAS 'Yeast' must be present.
    Only sucrose and glucose are taken into account for conversion.
    Conversion is based on reaction time, `tau`. Cleaning and unloading time,
    `tau_0`, fraction of working volume, `V_wf`, and number of reactors,
    `N_reactors`, are attributes that can be changed. Cost of a reactor
    is based on the NREL batch fermentation tank cost assuming volumetric
    scaling with a 6/10th exponent [1]_. 
    
    Parameters
    ----------
    ins : 
        Inlet fluids to be mixed into the fermentor.
    outs : 
        * [0] Vent
        * [1] Effluent
    tau : float
        Reaction time.
    N : int, optional
        Number of batch reactors
    V : float, optional
        Target volume of reactors [m^3].
    T=305.15 : float
        Temperature of reactor [K].
    P=101325 : float
        Operating pressure of reactor [Pa].
    Nmin=2 : int
        Minimum number of fermentors.
    Nmax=36: int
        Maximum number of fermentors.  
    efficiency=0.9 : float, optional
        User enforced efficiency.
    iskinetic=False: bool, optional
        If True, `Fermenation.kinetic_model` will be used.
    
    Notes
    -----
    Either N or V must be given.
    
    Examples
    --------
    Simulate a Fermentation object which models batch fermentation for the
    production of 1st generation ethanol using yeast.
    
    >>> from biorefineries.cane import create_sugarcane_chemicals
    >>> from biosteam.units import Fermentation
    >>> from biosteam import Stream, settings
    >>> settings.set_thermo(create_sugarcane_chemicals())
    >>> feed = Stream('feed',
    ...               Water=1.20e+05,
    ...               Glucose=1.89e+03,
    ...               Sucrose=2.14e+04,
    ...               DryYeast=1.03e+04,
    ...               units='kg/hr',
    ...               T=32+273.15)
    >>> F1 = NRELFermentation('F1',
    ...                   ins=feed, outs=('CO2', 'product'),
    ...                   tau=8, efficiency=0.90, N=8)
    >>> F1.simulate()
    >>> F1.show()
    NRELFermentation: F1
    ins...
    [0] feed
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.66e+03
                        Glucose  10.5
                        Sucrose  62.5
                        Yeast    456
    outs...
    [0] CO2
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    9.95
                        Ethanol  3.71
                        CO2      244
    [1] product
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.59e+03
                        Ethanol  240
                        Glucose  4.07
                        Yeast    532
    >>> F1.results()
    Fermentation                                       Units        F1
    Electricity         Power                             kW      66.6
                        Cost                          USD/hr       5.2
    Chilled water       Duty                           kJ/hr -1.41e+07
                        Flow                         kmol/hr  9.42e+03
                        Cost                          USD/hr      70.3
    Design              Reactor volume                    m3       247
                        Batch time                        hr      12.6
                        Loading time                      hr      1.57
                        Number of reactors                           8
                        Recirculation flow rate        m3/hr      17.7
                        Reactor duty                   kJ/hr -1.41e+07
                        Cleaning and unloading time       hr         3
                        Working volume fraction                    0.9
    Purchase cost       Heat exchangers (x8)             USD  1.57e+05
                        Reactors (x8)                    USD  1.87e+06
                        Agitators (x8)                   USD  1.17e+05
                        Cleaning in place                USD  8.89e+04
                        Recirculation pumps (x8)         USD  1.26e+05
    Total purchase cost                                  USD  2.36e+06
    Utility cost                                      USD/hr      75.5
    
    
    """
    line = 'Fermentation'
    _ins_size_is_fixed = False
    
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
    
    def _init(self, tau, N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
              efficiency=None, iskinetic=False, fermentation_reaction=None,
              cell_growth_reaction=None):
        NRELBatchBioreactor._init(self, tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self._load_components()
        self.iskinetic = iskinetic
        chemicals = self.chemicals
        self.hydrolysis_reaction = Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        if fermentation_reaction is None:
            if efficiency is None: efficiency = 0.9
            fermentation_reaction = Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', efficiency, chemicals)
        else:
            efficiency = fermentation_reaction.X
        self.fermentation_reaction = fermentation_reaction 
        if cell_growth_reaction is None:
            cell_growth_reaction = Reaction('Glucose -> Yeast', 'Glucose', 0.70, chemicals, basis='wt')
            cell_growth_reaction.basis = 'mol'
        self.cell_growth_reaction = cell_growth_reaction
        if all([i in self.chemicals for i in ('FFA', 'DAG', 'TAG', 'Glycerol')]):
            self.lipid_reaction = self.oil_reaction = ParallelReaction([
                Reaction('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
                Reaction('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
            ])
        else:
            self.lipid_reaction = None
        self.efficiency = efficiency
        
    def _calc_efficiency(self, feed, tau): # pragma: no cover
        # Get initial concentrations
        IDs = 'Yeast', 'Ethanol', 'Glucose', 
        X0, P0, S0 = feed.imass[IDs] / feed.F_vol
        
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
    def kinetic_model(z, t, *kinetic_constants): # pragma: no cover
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
        if P > Pm1: P = Pm1
        mu_X = mu_m1 * (S/(Ks1 + S)) * (1 - P/Pm1)**a*((1-X/Xm))
        mu_P = mu_m2 * (S/(Ks2 + S)) * (1 - P/Pm2)
        mu_S = mu_P / 0.45
        
        # Compute derivatives
        dXdt = mu_X * X
        dPdt = (mu_P * X)
        dSdt =  - mu_S * X
        return (dXdt, dPdt, dSdt)

    @property
    def efficiency(self):
        return self.fermentation_reaction.X
    @efficiency.setter
    def efficiency(self, efficiency):
        self.fermentation_reaction.X = efficiency

    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.hydrolysis_reaction.force_reaction(effluent)
        if self.iskinetic:
            self.fermentation_reaction.X = self._calc_efficiency(effluent, self._tau)
        self.fermentation_reaction.force_reaction(effluent)
        self.cell_growth_reaction.force_reaction(effluent)
        if self.lipid_reaction: self.lipid_reaction.force_reaction(effluent)
        effluent.empty_negative_flows()
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

BatchBioreactor = NRELBatchBioreactor
Fermentation = NRELFermentation