# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:57:12 2019

@author: yoelr
"""
from .._unit import Unit
from .decorators import cost

@cost('Flow rate', units='kg/hr', cost=2.5e6,
      CE=567.3, n=0.6, S=500e3, kW=3000, BM=1.39)
class Shredder(Unit):  pass