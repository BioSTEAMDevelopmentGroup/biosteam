# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:49 2019

@author: yoelr
"""
from .decorators import cost
from .._unit import Unit

__all__ = ('ConveyingBelt',)

@cost('Flow rate', CE=567, cost=813, ub=2000, n=0.38, N='Number of conveyors')
class ConveyingBelt(Unit):
    length = 40 #: ft
    height = 20 #: ft
    _N_outs = 1
    _has_power_utility = True
    _minimum_flow = 120
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        self.design_results['Flow rate'] = F_vol = feed.F_vol*35.315 # ft3/hr
        if F_vol < self._minimum_flow:
            self._lb_warning('Flow rate', F_vol, self._minimum_flow)
        F_mass = feed.F_mass*0.0006124 #lb/s
        self.power_utility(0.00058*F_mass**0.82*self.length
                           + self.height*0.00182*F_mass * 0.7457) # kW
        
        

