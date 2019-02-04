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
    """Create a Fermentation object which models batch fermentation for the production of 1st generation ethanol using yeast [1, 2]. Number of batches are selected based a loading time per volume ratio, *f*, assuming no influent surge time. Similation is based on reaction time, *tau*. Cleaning and unloading time, *tau_0*, and fraction of working volume, *V_wf*, are muttable properties. Cost of a reactor is based on the NREL batch fermentation tank cost assuming volumetric scaling with a 6/10th exponent [3].
    
    **Parameters**
    
        **tau:** Reaction time.
        
        **efficiency:** User enforced efficiency.
    
    **ins**
    
        [:] Inffluent streams
        
    **outs**
    
        [0] Effluent
        
        [1] CO2
    
    **References**
    
        [1] Samuel C. Oliveira, Dile P. Stremel, Eduardo C. Dechechi and Félix M. Pereira. Kinetic Modeling of 1‐G Ethanol Fermentations. Feb 2017. DOI: 10.5772/65460
        
        [2] S.C. Oliveira, H.F. De Castro, A.E.S. Visconti, R. Giudici. Continuous ethanol fermentation in a tower reactor with flocculating yeast recycle: scale-up effects on process performance, kinetic parameters and model predictions. Bioprocess Engineering 20 (1999) 525-530
        
        [3] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden National. Renewable Energy Laboratory Golden, Colorado. P. Schoen, J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical Report NREL/TP-5100-47764"""
        
    kwargs = {'tau': None, # Reaction time
              'efficiency': None}  # Theoretical efficiency
    _N_heat_util = 0
    _power_util = True
    
    #: Cleaning and unloading time
    tau_0 = 12 
    
    #: Loading time per volume (hr/m3)
    f = 0.001294
    
    #: Fraction of filled tank to total tank volume
    V_wf = 0.9 
    
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
    
    #: Base agitator power (hp)
    A_p = 30 
    
    #: tuple of kinetic parameters for the kinetic model. Default constants are fitted for Oliveria's model (mu_m1, mu_m2, Ks1, Ks2, Pm1, Pm2, Xm, Y_PS, a)
    kinetic_constants = (0.39,  # mu_m1
                         1.16,  # mu_m2
                         1.94,  # Ks1
                         2.6,   # Ks2
                         88.1,  # Pm1
                         126.5, # Pm2
                         114,   # Xm
                         0.45,  # Y_PS
                         0.16)  # a

    def _init(self):
        self._cooler = HXutility(ID=self.ID+' cooler', outs=())
    
    def _setup(self):
        if not self.kwargs['tau']:
            raise ValueError(f"Reaction time must be larger than 0, not '{self.kwargs['tau']}'")

    def _calc_efficiency(self, feed, tau):
        # Get initial concentrations
        y, e, s, w = feed.get_index(('Yeast',
                                     '64-17-5',
                                     '492-61-5',
                                     '7732-18-5'), CAS=True)
        mass = feed.mass
        volnet = feed.volnet
        concentration_in = mass/volnet
        X0, P0, S0 = [concentration_in[i] for i in (y, e, s)]
        X0 = 5
        
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
        eff = (S0 - Sf)/S0
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
        # 90% efficiency of fermentation
        # assume no Glycerol produced
        # no ethanol lost with CO2
        mw_Ethanol = 46.06844
        mw_glucose = 180.156

        # Unpack
        kwargs = self.kwargs
        out, CO2 = self.outs
        out.sum(out, self.ins)
        
        # Species involved in reaction
        species = ['Ethanol', 'Glucose', 'Sucrose', 'Water']
        e, g, s, w = indx = out.get_index(species)
        glucose = out.mol[indx[1]]
        sucrose = out.mol[indx[2]]
        
        # Hydrolysis
        new_glucose = sucrose*2
        ch_sucrose = -sucrose
        ch_Water = -sucrose
        out._mol[[g, s, w]] += [new_glucose, ch_sucrose, ch_Water]
        out.T = 32 + 273.15
        
        # Fermentation
        self.results['Operation']['Efficiency'] = eff = (
            kwargs.get('efficiency')
            or self._calc_efficiency(out, kwargs['tau'])
            )
        fermented = eff*(glucose + new_glucose)
        ch_glucose = - fermented
        ch_Ethanol = 0.5111*fermented * mw_glucose/mw_Ethanol
        co2 = out.get_index('CO2')
        changes = [ch_Ethanol, ch_glucose]
        out.mol[[e, g]] += changes
        CO2.mol[co2] = (1-0.5111)*fermented * mw_glucose/44.0095
        CO2.phase = 'g'
        CO2.T = out.T
    
    def _design(self):
        """
        * 'Number of reactors': (#)
        * 'Reactor volume':  (m3)
        * 'Cycle time': (hr)
        * 'Loading time': (hr)
        * 'Total dead time': (hr)
        """
        v_0 = self.outs[0].volnet
        f = self.f
        tau = self.kwargs['tau']
        tau_0 = self.tau_0
        Design = self.results['Design']
        Design.update(self._solve(v_0, tau, tau_0, f, self.V_wf))
        N = Design['Number of reactors']
        hx = self._cooler
        self.heat_utilities = hx.heat_utilities
        new_flow = Stream.like(self.outs[0])
        new_flow.mol /= N 
        hx.outs[0] = new_flow
        hx.ins[0] = new_flow
        self.heat_utilities[0](self.Hnet/N, self.outs[0].T)
        results = hx.heat_utilities[0].results
        hx._Duty = results['Duty']
        results['Duty'] *= N
        results['Cost'] *= N
        results['Flow'] *= N
        hx._design()
        hx._cost()
        return Design
    
    def _cost(self):
        """
        * 'Reactors': (USD)
        * 'Coolers': (USD)
        * 'Agitators': (USD)
        """
        results = self.results
        Cost = results['Cost']
        Design = results['Design']
        N = Design['Number of reactors']
        V_i = Design['Reactor volume']
        R = V_i/self.V_0
        F_CEPCI = self.CEPCI/self.CEPCI_0
        FF = F_CEPCI * N * R**self.exp
        
        Cost['Reactors'] = FF * self.C_0
        Cost['Coolers'] = N * self._cooler.results['Cost']['Heat exchanger']
        self.power_utility(FF * self.A_p)
        Cost['Agitators'] = FF * self.C_A
        return Cost
