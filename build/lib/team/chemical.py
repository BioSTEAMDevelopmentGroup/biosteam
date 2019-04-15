# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:50:03 2018

@author: yoelr
"""
from thermo.chemical import Chemical as TChem
from .compound import Compound

__all__ = ('Chemical', 'Solid', 'Liquid', 'Gas')

class Chemical(Compound, TChem):
    """An extension of the ChEDL thermo Chemical class. The enthalpy property, 'H', does not account for excess ethalpy, it is simply based on latent heats and heat capacities at a constant pressure of 101325. All thermodynamic energies are now calculated properties and not attributes. 

`Read the docs for thermo.Chemical for accurate documentation. <http://thermo.readthedocs.io/en/latest/thermo.chemical.html>`__"""

    def __init__(self, ID, T=298.15, P=101325):
        ID = ID.replace('_', ' ')
        TChem.__init__(self, ID, T, P)
        self.ID = ID.replace(' ', '_')
        if self.CAS == '56-81-5':
            self.__UNIFAC_Dortmund_groups = {2: 2, 3: 1, 14: 2, 81: 1}
        self.Hfm = self.Hf

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
        try:
            self._eos_T_101325 = self.eos.to_TP(self.T, 101325)
        except:
            pass


class Solid(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as solid."""
    
    def __init__(self, ID, T=298.15, P=101325):
        super().__init__(ID, T, P)
    
    @property
    def phase(self):
        return 's'
    @phase.setter
    def phase(self, phase): pass


class Liquid(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as liquid."""
    
    def __init__(self, ID, T=298.15, P=101325):
        super().__init__(ID, T, P)
    
    @property
    def phase(self):
        return 'l'
    @phase.setter
    def phase(self, phase): pass


class Gas(Chemical):
    """Create a :doc:`Chemical` such that its phase remains as gas."""
    
    def __init__(self, ID, T=298.15, P=101325):
        super().__init__(ID, T, P)
    
    @property
    def phase(self):
        return 'g'
    @phase.setter
    def phase(self, phase): pass





