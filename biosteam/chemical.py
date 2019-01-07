# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:50:03 2018

@author: yoelr
"""
from thermo.chemical import Chemical as TChem
from .specie import Specie


class Chemical(Specie, TChem):
    """An extension on the ChEDL thermo Chemical class through inheritance of the Specie class. The enthalpy property, 'H', now ignores excess ethalpies and calculates based on latent heats and heat capacities at a constant pressure of 101325. All thermodynamic energies are now calculated properties and not attributes.

`Read the docs for thermo.Chemical for accurate documentation. <http://thermo.readthedocs.io/en/latest/thermo.chemical.html>`__"""

    def __init__(self, ID, T=298.15, P=101325):
        TChem.__init__(self, ID, T, P)
        self.ID = ID
        if self.CAS == '56-81-5':
            self.__UNIFAC_Dortmund_groups = {2: 2, 3: 1, 14: 2, 81: 1}

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