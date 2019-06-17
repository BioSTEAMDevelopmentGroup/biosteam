# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:40:50 2019

@author: yoelr
"""
from .. import Unit
from .decorators import cost
from .metaclasses import static

@cost('Flow rate', limit=10e4, CE=567, cost=1096, exp=0.22)
class ScrewFeeder(Unit, metaclass=static):
    length = 30 #: ft
    _N_outs = 1
    _has_power_utility = True
    _minimum_flow = 400
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        r = self.results
        Design = r['Design']
        volnet = feed.volnet*35.315 # ft3/hr
        if volnet < self._minimum_flow:
            self._lb_warning('Flow rate', volnet, self._minimum_flow)
        Design['Flow rate'] = volnet
        massnet = feed.massnet*0.0006124 #lb/s
        self._power_utility(0.0146*massnet**0.85*self.length*0.7457)