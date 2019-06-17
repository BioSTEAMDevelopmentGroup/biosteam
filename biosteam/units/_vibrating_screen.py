# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 20:51:22 2019

@author: yoelr
"""
from .. import Unit
from .decorators import design, cost
from .metaclasses import splitter

@cost('Area', limit=200, CE=567, cost=1010, exp=0.91)
@design('Area', 'ft^2', lambda self: self.ins[0].massnet/(self.capacity*self.mesh_opening))
class VibratingScreen(Unit, metaclass=splitter):
    # Assume 3-deck vibrating screen
    
    #: Flow rate per area of screen per apareture (kg∕ft2-hr-mm)
    capacity = 0.2*907.18474
    
    #: Mesh opening (mm)
    mesh_opening = 8
    