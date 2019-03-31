# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:45:47 2018

@author: yoelr
"""
from biosteam.units.reactor import BatchReactor
from biosteam.units.hx import HXutility
from biosteam import Stream, np
from scipy.integrate import odeint

class Fermentation(BatchReactor):
    """
    Create a Fermentation object which models large-scale batch fermentation for the production of 1st generation ethanol using yeast [1, 2, 3, 4]. Only sucrose and glucose are taken into account. Conversion is based on reaction time, *tau*. Cleaning and unloading time, *tau_0*, fraction of working volume, *V_wf*, and number of reactors, N_reactors, are attributes that can be changed. Cost of a reactor is based on the NREL batch fermentation tank cost assuming volumetric scaling with a 6/10th exponent [3].
    
    **Parameters**
    
        **tau:** Reaction time.
        
        **efficiency:** User enforced efficiency.
    
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
    _kwargs = {'tau': None, # Reaction time
               'efficiency': None}  # Theoretical efficiency
    _N_heat_utilities = 0
    _has_power_utility = True
    
    #: Cleaning and unloading time
    tau_0 = 12 
    
    #: Number of reactors
    _N_reactors = None
    
    #: Fraction of filled tank to total tank volume
    V_wf = 0.8
    
    #: Base CEPCI 
    CEPCI_0 = 521.9 
    
    #: Base price (USD)    
    C_0 = 844000
    
    #: Base volume (m3)    
    V_0 = 3785
    
    #: Scaling exponent for costing
    exp = 0.6 
    
    #: Base agitator price (USD)
    C_A = 52500 
    
    #: Base agitator power (kW)
    A_p = 22.371
    
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
    @property
    def N_reactors(self):
        """Enforced Number of batch reactors."""
        return self._N_reactors
    
    @N_reactors.setter
    def N_reactors(self, N):
        if N <= 1: raise ValueError(f"Number of reactors must be greater than 1, value {N} is infeasible")
        self._N_reactors = N

    def _init(self):
        self._cooler = hx = HXutility('*')
        self.heat_utilities = hx.heat_utilities
    
    def _setup(self):
        if not self._kwargs['tau']:
            raise ValueError(f"Reaction time must be larger than 0, not '{self._kwargs['tau']}'")

    def _calc_efficiency(self, feed, tau):
        # Get initial concentrations
        y, e, s, w = feed.indices('Yeast',
                                  '64-17-5',
                                  '492-61-5',
                                  '7732-18-5', CAS=True)
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
        

    def _run(self):
        # Assume no Glycerol produced
        # No ethanol lost with CO2
        mw_Ethanol = 46.06844
        mw_glucose = 180.156

        # Unpack
        kwargs = self._kwargs
        out, CO2 = self._outs
        out.sum(out, self._ins)
        
        # Species involved in reaction
        species = ['64-17-5', '492-61-5', '57-50-1', '7732-18-5']
        e, g, s, w = indx = out.indices(*species, CAS=True)
        glucose = out.mol[indx[1]]
        sucrose = out.mol[indx[2]]
        
        # Hydrolysis
        new_glucose = sucrose*2
        ch_sucrose = -sucrose
        ch_Water = -sucrose
        out._mol[[g, s, w]] += [new_glucose, ch_sucrose, ch_Water]
        out.T = 32 + 273.15
        
        # Fermentation
        self._results['Design']['Efficiency'] = eff = (
            kwargs.get('efficiency')
            or self._calc_efficiency(out, kwargs['tau'])
            )
        fermented = eff*(glucose + new_glucose)
        ch_glucose = - fermented
        ch_Ethanol = 0.5111*fermented * mw_glucose/mw_Ethanol
        co2 = out.indices('124-38-9', CAS=True)
        changes = [ch_Ethanol, ch_glucose]
        out.mol[[e, g]] += changes
        CO2.mol[co2] = (1-0.5111)*fermented * mw_glucose/44.0095
        CO2.phase = 'g'
        CO2.T = out.T
    
    def _optimize_N_reactors(self):
        Summary = self._results['Summary']
        cost_old = np.inf
        self.N_reactors = 2
        self.simulate()
        cost_new = Summary['Purchase cost']
        while cost_new < cost_old:
            self.N_reactors += 1
            self.simulate()
            cost_old = cost_new
            cost_new = Summary['Purchase cost']
        self.N_reactors -= 1
        
    def _design(self):
        N_reactors = self.N_reactors
        if N_reactors is None:
            self._N_reactors = 2
            self._optimize_N_reactors()
            N_reactors = self._N_reactors
            self._N_reactors = None
            
        v_0 = self.outs[0].volnet
        tau = self._kwargs['tau']
        tau_0 = self.tau_0
        Design = self._results['Design']
        Design['Number of reactors'] = N_reactors
        Design.update(self._solve(v_0, tau, tau_0, N_reactors, self.V_wf))
        hx = self._cooler
        new_flow = Stream.like(self.outs[0])
        new_flow.mol /= N_reactors 
        hx.outs[0] = new_flow
        hx.ins[0] = new_flow
        hu = hx.heat_utilities[0]
        hu(self.Hnet/N_reactors, self.outs[0].T)
        hx._design(hu.duty)
        hx._cost()
        hu.duty *= N_reactors
        hu.cost *= N_reactors
        hu.flow *= N_reactors
        return Design
    
    def _cost(self):
        results = self._results
        Cost = results['Cost']
        Design = results['Design']
        N = Design['Number of reactors']
        V_i = Design['Reactor volume']
        R = V_i/self.V_0
        F_CEPCI = self.CEPCI/self.CEPCI_0
        FF = F_CEPCI * N * R**self.exp
        
        Cost['Reactors'] = FF * self.C_0
        Cost['Coolers'] = N * self._cooler._results['Cost']['Heat exchanger']
        self.power_utility(self.A_p * R * N * F_CEPCI)
        Cost['Agitators'] = FF * self.C_A
        return Cost
