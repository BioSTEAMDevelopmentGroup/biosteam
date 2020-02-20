# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:32:01 2019

@author: yoelr
"""

from .._unit import Unit
from .decorators import cost

__all__ = ('MagneticSeparator',)

@cost('Flow rate', units='kg/hr', CE=576, cost=533471, S=333333, n=0.6)
class MagneticSeparator(Unit): pass
    