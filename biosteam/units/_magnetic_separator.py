# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:32:01 2019

@author: yoelr
"""

from .. import Unit
from .decorators import cost, design

@cost('Flow rate', CE=576, cost=533471, S=333333, exp=0.6)
@design('Flow rate', 'kg/hr', lambda self: self.ins[0].massnet)
class MagneticSeparator(Unit):
    _N_outs = 1
    _linkedstreams = True
    