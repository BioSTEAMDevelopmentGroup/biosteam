# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .decorators import cost
from .. import Unit

__all__ = ('VentScrubber',)

@cost('Flow rate', units='kg/hr',
      S=22608, CE=522, cost=215e3, n=0.6, BM=2.4)
class VentScrubber(Unit): 
    _N_ins = _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, gas):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.gas = gas
    
    def _run(self):
        water, vent_entry = self.ins
        vent_exit, bottoms = self.outs
        vent_exit.copy_like(vent_entry)
        bottoms.empty()
        bottoms.copy_flow(vent_exit, self.gas,
                          remove=True, exclude=True)
        bottoms.mix_from([bottoms, water], energy_balance=False)
        
    def _design(self):
        self.design_results['Flow rate'] = self._ins[1].F_mass