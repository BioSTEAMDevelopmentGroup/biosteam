# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:57:12 2019

@author: yoelr
"""
from .. import Unit
from .decorators import cost, design
from .metaclasses import static

@cost('Flow rate', cost=2.5e6, CE=567.3, exp=0.6, S=500e3, kW=3000)
@design('Flow rate', 'kg/hr')
class Shredder(Unit, metaclass=static): pass