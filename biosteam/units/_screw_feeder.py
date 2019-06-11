# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:40:50 2019

@author: yoelr
"""
from .. import Unit
import numpy as np
from .decorators import cost

@cost('Flow rate', N='N', CE=567, cost=1096, exp=0.22)
class ScrewFeeder(Unit):
    length = 30 #: ft
    _N_outs = 1
    _has_power_utility = True
    _linkedstreams = True
    _bounds = {'Volumetric flow': (400, 10000)}
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        r = self.results
        Design = r['Design']
        volbounds = self.bounds['Volumetric flow']
        volnet = feed.volnet*35.315 # ft3/hr
        N = volnet/volbounds[-1]
        Design['N'] = N = np.ceil(N)
        Design['Flow rate'] = volnet
        massnet = feed.massnet*0.0006124 #lb/s
        self._power_utility(0.0146*massnet**0.85*self.length*0.7457)