# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 19:27:01 2018

@author: yoelr
"""
from ._compound import Compound
import numpy as np

__all__ = ('Substance',)


def select_value(obj, user_defined, attr, default):
    if user_defined is not None: return user_defined
    try:
        obj_value = getattr(obj, attr)
        if obj_value is None:
            return default
        else:
            return obj_value
    except:
        return default
    
# %%

class Substance(Compound):
    """Create a Substance object with the same material properties as "obj" argument for simplified modeling purposes. The phase stays constant as the specified `default_phase`.

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
        
        **default_phase:** {'l', 's', 'g'} Phase of substance
        
    **Examples**
    
        Create a 'Yeast' substance with the same properties as water:
    
        .. code-block:: python
            
            >>> from biosteam.compounds import Substance
            >>> Yeast = Substance('Yeast', obj=Chemical('Water'))
            >>> Yeast.rho
            997
            
    """
    def __init__(self, ID, CAS='', obj=None, MW=1, T=298.15, P=101325,
                 Cp=None, rho=None, k=None, mu=None, sigma=None,
                 Pr=None, Hfm=None, Hc=None, obj_state = None,
                 default_phase='l'):
        if obj:
            obj.T = T
            obj.P = P
            obj.phase = default_phase
            if obj_state:
                for i, j in obj_state.items():
                    setattr(obj, i, j)
        
        self.ID = ID
        self.T = T
        self.P = P
        self.MW = MW
        self.Hvapm = 0
        
        # Set property values (defaults to Water property values)
        self.CAS = select_value(obj, CAS, 'CAS', ID)
        self.Cp = select_value(obj, Cp, 'Cp', 4.18)
        self.sigma = select_value(obj, sigma, 'mu', 0.072055)
        self.rho = select_value(obj, rho, 'rho', 997)
        self.mu = select_value(obj, mu, 'mu', 0.00091272)
        self.k = select_value(obj, k, 'k', 0.5942)
        self.Hfm = select_value(obj, Hfm, 'Hfm', 0)
        self.Hc = select_value(obj, Hc, 'Hc', 0)
        self.rhom = self.rho/MW*1000
        self.Vm = 1/self.rhom
        self.Cpm = self.Cp * MW
        self.alpha = self.k/self.rho/self.Cp
        self.nu = self.mu/self.rho
        self.Pr = self.nu/self.alpha
        self.UNIFAC_groups = self.UNIFAC_Dortmund_groups = None
        self.phase_ref = default_phase
        
        if default_phase not in 'slg':
            raise ValueError(f"default phase must be either 's', 'l' or 'g', not '{default_phase}'")
        
        # Phase change temperature consitent with phase
        self.Tb = 0 if default_phase=='g' else np.inf 
        self.Tm = np.inf if default_phase=='s' else 0
        
    @property
    def Cplm(self):
        return self.Cpm
    @property
    def phase(self):
        """Phase of substance."""
        return self.phase_ref
    @phase.setter
    def phase(self, phase): pass