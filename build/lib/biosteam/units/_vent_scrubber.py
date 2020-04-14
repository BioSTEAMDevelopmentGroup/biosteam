# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 01:11:46 2019

@author: yoelr
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
        bottoms.copy_flow(vent_exit, self.gas,
                         remove=True, exclude=True)
        bottoms.mol[:] += water.mol
        
    def _design(self):
        self.design_results['Flow rate'] = self._outs[0].F_mass