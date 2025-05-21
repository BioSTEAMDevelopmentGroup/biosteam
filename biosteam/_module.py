# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from ._facility import Facility

__all__ = ('Module', 'FacilityModule')

class Module(Unit):
    _N_ins = 0
    _N_outs = 0
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def _assert_compatible_property_package(self): pass
    
    def _setup(self):
        self.auxiliary_system._setup(load_configuration=False)
        super()._setup()
    
    def _run(self):
        self.auxiliary_system.converge()
    
    def _design(self):
        for i in self.auxiliary_system.units: i._design()

    def _cost(self):
        for i in self.auxiliary_system.units: i._cost()
        
class FacilityModule(Module, Facility): network_priority = None