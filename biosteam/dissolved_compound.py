# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 19:27:01 2018

@author: yoelr
"""
from biosteam.compound import Compound
from biosteam import np

class DissolvedCompound(Compound):
    """
    Creates a DissolvedCompound object with the same material properties as "obj" argument for simplified modeling purposes. The phase stays constant as a liquid ('l').

    **Parameters**
    
        **ID:** [str] ID of specie
        
        **CAS:** [str] CAS identification
        
        **obj** [Specie] Should contain analog properties for specie at STP. Defaults to water if None.

        **MW:** Molecular weight of specie
         
        **T:** Temperature (K)
         
        **P:** Pressure (Pa)
         
        **rho:** Density (kg/m^3)
        
        **k:** Thermal conductivity (W/m/K)
        
        **mu:** Hydrolic viscosity (Pa*s)
        
        **sigma:** Surface tension (N/m)
        
        **Hfm:** Heat of formation (J/mol)
        
        **Hc:** Heat of combustion (J/mol)
        
    **Examples**
    
        Create a 'Yeast' compound with the same properties as water:
    
        .. code-block:: python
            
            >>> from biosteam import DissolvedCompound
            >>> Yeast = DissolvedCompound('Yeast', obj=Chemical('Water'))
            >>> Yeast.rho
            997
            
    """
    Tb = np.inf #: No boiling point (K)
    phase_ref = 'l'
    
    def __init__(self, ID, CAS='', obj=None, MW=1, T=298.15, P=101325,
                 Cp=None, rho=None, k=None, mu=None,
                 sigma=None, Pr=None, Hfm=None, Hc=None):
        self.ID = ID
        self.T = T
        self.P = P
        self.MW = MW
        self.Hvapm = 0
        
        # Default values
        Cp_d = 4.18
        rho_d = 997
        k_d = 0.5942
        mu_d = 0.00091272
        sigma_d = 0.072055
        Hfm_d = 0
        Hc_d = 0
        
        # Set object default values
        if obj:
            Cp_d = obj.Cp if obj.Cp else Cp_d
            rho_d = obj.rho if obj.rho else rho_d
            k_d = obj.k if obj.k else k_d
            mu_d = obj.mu if obj.mu else mu_d
            sigma_d = obj.sigma if obj.sigma else sigma_d
            Hfm_d = obj.Hf if obj.Hf else Hfm_d
            Hc_d = obj.Hc if obj.Hc else Hc_d
        
        # Set all values
        self.CAS = CAS if CAS else ID
        self.Cp = Cp if Cp else Cp_d
        self.Cpm = self.Cp * MW
        self.rho = rho if rho else rho_d
        self.rhom = self.rho/MW*1000
        self.Vm = 1/self.rhom
        self.k = k if k else k_d
        self.alpha = self.k/self.rho/self.Cp
        self.mu = mu if mu else mu_d
        self.nu = self.mu/self.rho
        self.sigma = sigma if sigma else sigma_d
        self.Pr = self.nu/self.alpha
        self.UNIFAC_groups = {}
        self.Hfm = Hfm if Hfm else Hfm_d
        self.Hc = Hc if Hc else Hc_d
    
    @property
    def Cplm(self):
        return self.Cpm
    @property
    def phase(self):
        """Phase is always 'l' (liquid), regardless."""
        return 'l'
    @phase.setter
    def phase(self, phase): pass