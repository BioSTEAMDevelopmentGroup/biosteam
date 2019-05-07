# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:57:12 2019

@author: yoelr
"""
from .. import Unit
from .decorators import cost, design

@cost('Flow rate', cost=2.5e6, CE=567.3, exp=0.6, S=500e3, kW=0.006)
@design('Flow rate', 'kg/hr', lambda self: self._ins[0].massnet)
class Shredder(Unit):
    _linkedstreams = True
    _N_outs = 1