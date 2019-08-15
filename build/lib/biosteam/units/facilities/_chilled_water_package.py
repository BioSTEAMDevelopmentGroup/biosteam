# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ..decorators import cost

__all__ = ('ChilledWaterPackage',)

@cost('Duty', S=-14*4184000, kW=3400*0.7457, cost=1375e3, CE=551, n=0.7, BM=1.5)
class ChilledWaterPackage(Facility):
    """Create a chilled water package that is cost based on flow rate of chilled water.
    
    **Parameters**
    
        **Duty:** Chilled water duty (kJ/hr)
    
    """
    _N_heat_utilities = 1
    _N_ins = _N_outs = 0
    _units = {'Duty': 'kJ/hr'}
    def __init__(self, ID=''):
        Facility.__init__(self, ID, None, None)
        self.chilled_water_utilities = set()
        
    def _design(self):
        cwu = self.chilled_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u._heat_utilities:
                    if hu.ID == 'Chilled water': cwu.add(hu)
        self._Design['Duty'] = duty = sum([i.duty for i in cwu])
        self._heat_utilities[0](duty, 330)
        
        