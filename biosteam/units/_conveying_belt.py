# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:49 2019

@author: yoelr
"""
from .decorators import cost
from ._static import Static

@cost('Flow rate', CE=567, cost=813, ub=500, n=0.38)
class ConveyingBelt(Static):
    length = 40 #: ft
    height = 20 #: ft
    _N_outs = 1
    _has_power_utility = True
    _minimum_flow = 120
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        self._Design['Flow rate'] = volnet = feed.volnet*35.315 # ft3/hr
        if volnet < self._minimum_flow:
            self._lb_warning('Flow rate', volnet, self._minimum_flow)
        massnet = feed.massnet*0.0006124 #lb/s
        self._power_utility(0.00058*massnet**0.82*self.length + self.height*0.00182*massnet * 0.7457) # kW
        
        

