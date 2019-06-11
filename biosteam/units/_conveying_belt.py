# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:10:49 2019

@author: yoelr
"""
from .. import Unit
from .decorators import cost
import numpy as np

@cost('Flow rate', N='N', CE=567, cost=813, exp=0.38)
class ConveyingBelt(Unit):
    length = 40 #: ft
    height = 20 #: ft
    _N_outs = 1
    _has_power_utility = True
    _linkedstreams = True
    _bounds = {'Flow rate': (120, 500)}
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        Design = self._results['Design']
        volbounds = self._bounds['Flow rate']
        volnet = feed.volnet*35.315 # ft3/hr
        Design['N'] = N = np.ceil(volnet/volbounds[-1])
        Design['Flow rate'] = volnet/N
        massnet = feed.massnet*0.0006124 #lb/s
        self._power_utility(0.00058*massnet**0.82*self.length + self.height*0.00182*massnet * 0.7457) # kW
        
        

