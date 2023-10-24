# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
from ._batch_bioreactor import BatchBioreactor
from scipy.integrate import odeint
from thermosteam.reaction import Reaction, ParallelReaction

__all__ = ('Fermentation',)

class Fermentation(BatchBioreactor):
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
    >>> F1 = Fermentation('F1',
    ...                   ins=feed, outs=('CO2', 'product'),
    ...                   tau=8, efficiency=0.90, N=8)
    >>> F1.simulate()
    >>> F1.show()
    Fermentation: F1
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
    
    References
    ----------
    .. [1] Oliveira, Samuel C., et al. "Discrimination between ethanol 
        inhibition models in a continuous alcoholic fermentation process using
        flocculating yeast." Applied biochemistry and biotechnology 74.3 (1998): 161-172.
    
    .. [2] Oliveira, Samuel C., et al. "Continuous ethanol fermentation in a
        tower reactor with flocculating yeast recycle: scale-up effects on process
        performance, kinetic parameters and model predictions." Bioprocess
        Engineering 20.6 (1999): 525-530.
    
    .. [3] Oliveira, Samuel C., et al. "Mathematical modeling of a continuous
        alcoholic fermentation process in a two-stage tower reactor cascade with
        flocculating yeast recycle." Bioprocess and biosystems engineering 38.3
        (2015): 469-479.
    
    .. [4] Oliveira, Samuel C., et al. "Kinetic Modeling of 1‐G Ethanol
        Fermentations." Fermentation Processes. InTech, 2017.
    
    .. [5] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
        National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
        J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
        Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
        for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
        Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
        Report NREL/TP-5100-47764
    
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
        BatchBioreactor._init(self, tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
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
