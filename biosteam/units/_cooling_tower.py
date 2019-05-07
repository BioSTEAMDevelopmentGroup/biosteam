# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from biosteam import Unit
from .decorators import cost, design

__all__ = ('CoolingTower',)

@cost('Flow rate', cost=1100, CE=567, exp=0.68)
@design('Flow rate', 'gpm', lambda self: 0.07934*self.ins[0].molnet)
class CoolingTower(Unit):
    _N_ins = 1; _N_outs = 0