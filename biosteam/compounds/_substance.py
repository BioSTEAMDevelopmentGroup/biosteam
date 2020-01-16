# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 19:27:01 2018

@author: yoelr
"""
from ._compound import Compound, _R
from ._chemical import Chemical
import numpy as np
from math import log

__all__ = ('Substance', 'StaticChemical')


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

    Parameters
    ----------
    ID : str
         ID of specie
    CAS="" : str
        CAS identification
    obj=None : Comound, defaults to Chemical('Water')
        Should contain analog properties for substance at STP.
    MW=1 : float
        Molecular weight of compound.
    T=298.15 : float
        Temperature (K)
    P=101325 : float
        Pressure (Pa)
    rho=None : float, defaults to obj.rho
        Density (kg/m^3)
    k=None : float, defaults to obj.k
        Thermal conductivity (W/m/K)
    mu=None : float, defaults to obj.mu
        Hydrolic viscosity (Pa*s)
    sigma=None : float, defaults to obj.sigma
        Surface tension (N/m)
    Hfm : float, defaults to obj.Hfm
        Heat of formation (J/mol)
    Hc : float, defaults to obj.Hc
        Heat of combustion (J/mol)
    phase=None : {'l', 's', 'g'}, defaults to obj.phase
        Phase of substance.
    **properties
        Additional properties to set.
        
    Examples
    --------
    Create a 'Yeast' substance with the same properties as water:

    .. code-block:: python
        
        >>> from biosteam.compounds import Substance
        >>> Yeast = Substance('Yeast', obj=Chemical('Water'))
        >>> Yeast.rho
        997
            
    """
    def __init__(self, ID, CAS='', obj=None, *,
                 MW=1, T=298.15, P=101325, Cp=None,
                 rho=None, k=None, mu=None, sigma=None,
                 Pr=None, Hf=None, Hc=None, phase=None, **properties):
        if obj:
            obj.T = T
            obj.P = P
            if phase:
                obj.phase = phase
            else:
                phase = obj.phase
        if not phase: phase = 'l'
        
        self.ID = ID.replace(' ', '_')
        self.T = T
        self.P = P
        self.Hvapm = 0
        
        # Set property values (defaults to Water property values)
        self.MW = MW = select_value(obj, MW, 'MW', 1)
        self.CAS = select_value(obj, CAS, 'CAS', ID)
        self.Cp = select_value(obj, Cp, 'Cp', 4.18)
        self.sigma = select_value(obj, sigma, 'mu', 0.072055)
        self.rho = select_value(obj, rho, 'rho', 997)
        self.mu = select_value(obj, mu, 'mu', 0.00091272)
        self.k = select_value(obj, k, 'k', 0.5942)
        self.Hf = select_value(obj, Hf, 'Hf', 0)
        self.Hc = select_value(obj, Hc, 'Hc', 0)
        self.Hfm = self.Hf*MW
        self.Hcm = self.Hc*MW
        self.rhom = self.rho/MW*1000
        self.Vm = 1/self.rhom
        self.Cpm = self.Cp * MW
        self.alpha = self.k/self.rho/self.Cp
        self.nu = self.mu/self.rho
        self.Pr = self.nu/self.alpha
        self.UNIFAC_groups = self.UNIFAC_Dortmund_groups = None
        self.phase_ref = phase
        
        if phase not in 'slg':
            raise ValueError(f"default phase must be either 's', 'l' or 'g', not '{phase}'")
        
        # Phase change temperature consitent with phase
        self.Tb = 0 if phase=='g' else np.inf 
        self.Tm = np.inf if phase=='s' else 0
        for i in properties: setattr(self, i, properties[i])
        
    def copy(self):
        copy = self.__new__(self.__class__)
        copy.__dict__.update(self.__dict__)
        return copy
        
    __copy__ = copy
    
    @property
    def H(self):
        return self.Cpm*(self.T-self.T_ref)
    @property
    def S(self):
        return self.Cpm*log(self.T/self.T_ref) + (_R*log(self.P/self.P_ref) if self.phase=='g' else 0)
    
    @property
    def phase(self):
        """Phase of substance."""
        return self.phase_ref
    @phase.setter
    def phase(self, phase): pass


def StaticChemical(ID, *, T=298.15, P=101325, MW=None,
                   Cp=None, rho=None, k=None, mu=None, sigma=None,
                   Pr=None, Hf=None, Hc=None, phase=None, **properties):
    obj = Chemical(ID, T, P)
    return Substance(obj=obj, ID=obj.ID, CAS=obj.CAS,
                     MW=obj.MW, T=T, P=P, Cp=Cp, rho=rho,
                     k=k, mu=mu, sigma=sigma, Pr=Pr, Hf=Hf,
                     Hc=Hc, phase=phase, **properties)