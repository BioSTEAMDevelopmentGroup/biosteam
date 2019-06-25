# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:32:01 2019

@author: yoelr
"""

from .. import Unit
from .decorators import cost
from .metaclasses import static

@cost('Flow rate', units='kg/hr', CE=576, cost=533471, S=333333, n=0.6)
class MagneticSeparator(Unit, metaclass=static): pass
    