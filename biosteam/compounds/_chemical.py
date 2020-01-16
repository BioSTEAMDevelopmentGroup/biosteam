# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:50:03 2018

@author: yoelr
"""
from ..thermo.chemical import Chemical as TChem
from ._compound import Compound
import numpy as np
import re

__all__ = ('Chemical', 'Solid', 'Liquid', 'Gas')

class Chemical(Compound, TChem):
    """An extension of the ChEDL thermo Chemical class. The enthalpy property, 'H', does not account for excess ethalpy, it is simply based on latent heats and heat capacities at a constant pressure of 101325. All thermodynamic energies are now calculated properties and not attributes. 

`Read the docs for thermo.Chemical for accurate documentation. <http://thermo.readthedocs.io/en/latest/thermo.chemical.html>`__"""

    def __new__(cls, ID, T=298.15, P=101325):
        try:
            search_ID = re.sub(r"\B([A-Z])", r" \1", ID).capitalize().replace('_', ' ')
            self = TChem.__new__(cls, search_ID, T, P)
        except:
            raise LookupError(f"chemical '{ID}' not found in data bank, try chemical's CAS")
        self.ID = ID.replace(' ', '_')
        if self.CAS == '56-81-5':
            self.__UNIFAC_Dortmund_groups = {2: 2, 3: 1, 14: 2, 81: 1}
        self.Hfm = self.Hf
        return self

    def __copy__(self):
        copy = object.__new__(self.__class__)
        copy.__dict__.update(self.__dict__)
        return copy
    copy = __copy__

    # These properties are gone
    Hm = None
    Sm = None
    Um = None
    Am = None
    
    # New units for heat capacity
    @property
    def Cp(self):
        "Specific heat capacity (kJ/kg/K)"
        return self.Cpm/self.MW

    # Baseline equation of state set at atmospheric pressure
    def set_thermo(self):
        try: self._eos_T_101325 = self.eos.to_TP(self.T, 101325)
        except: pass


class Solid(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as solid."""
    def __new__(cls, ID, T=298.15, P=101325):
        self = super().__new__(cls, ID, T, P)
        # Phase change temperature consitent with phase
        self.Tb = np.inf 
        self.Tm = np.inf
        self.phase_ref = 's'
        return self
    
    @property
    def phase(self): return 's'
    @phase.setter
    def phase(self, phase): pass


class Liquid(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as liquid."""
    phase_ref = 'l'
    def __new__(cls, ID, T=298.15, P=101325):
        self = super().__new__(cls, ID, T, P)
        # Phase change temperature consitent with phase
        self.Tb = np.inf 
        self.Tm = 0
        self.phase_ref = 'l'
        return self
    
    @property
    def phase(self): return 'l'
    @phase.setter
    def phase(self, phase): pass


class Gas(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as gas."""
    phase_ref = 'g'
    def __new__(cls, ID, T=298.15, P=101325):
        self = super().__new__(cls, ID, T, P)
        # Phase change temperature consitent with phase
        self.Tb = 0
        self.Tm = 0
        self.phase_ref = 'g'
        return self
    
    @property
    def phase(self): return 'g'
    @phase.setter
    def phase(self, phase): pass



