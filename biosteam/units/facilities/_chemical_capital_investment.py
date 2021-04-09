# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from . import Facility

__all__ = ('ChemicalCapitalInvestment',)

class ChemicalCapitalInvestment(Facility):
    """
    Create a ChemicalCapitalInvestment facility object that computes
    the volume of a chemical required in a process and the incurred capital cost.
    
    Parameters
    ----------
    ID : str
        ID of unit operation.
    chemical : str
        ID of chemical.
    price : flaot
        Price of chemical [USD / m^3].
    
    """
    _F_BM_default = {'Chemical': 1.}
    _N_heat_utilities = 0
    network_priority = 3
    _N_outs = _N_ins = 0
    
    def __init__(self, ID, chemical, price):
        self.chemical = chemical
        self.price = price
        super().__init__(ID)
        
    def _run(self): pass

    def _design(self):
        self.design_results['Volume'] = bst.process_tools.volume_of_chemical_in_units(self.system.units, self.chemical)
        
    def _cost(self):
        self.purchase_costs['Chemical'] = self.price * self.design_results['Volume']