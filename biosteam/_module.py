# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .utils import piping
from ._unit import Unit
from ._facility import Facility

__all__ = ('Module', 'FacilityModule')

class Module(Unit):
    _ins_size_is_fixed = _outs_size_is_fixed = False
    _N_ins = _N_outs = 0
    
    def _init(self, system, run_full_simulation=None):
        if not self._ins: self._ins = piping.StreamPorts.from_inlets(system.feeds)
        if not self._outs: self._outs = piping.StreamPorts.from_outlets(system.products)
        self.register_auxiliary(system, 'auxsystem')
        if run_full_simulation is None: run_full_simulation = bool(system.facilities)
        self.run_full_simulation = run_full_simulation
    
    def _assert_compatible_property_package(self): pass
    
    def _setup(self):
        self.auxsystem._setup(load_configuration=False)
        super()._setup()    
    
    def _run(self):
        if self.run_full_simulation:
            self.auxsystem.simulate()
        else:
            self.auxsystem.converge()

        
class FacilityModule(Module, Facility): network_priority = None