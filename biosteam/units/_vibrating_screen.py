# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 20:51:22 2019

@author: yoelr
"""
from .. import Unit
from ._splitter import Splitter
from numpy import ceil
from .decorators import cost

@cost('Area', 'Vibrating screens', N='N', CE=567, cost=1010, exp=0.91)
class VibratingScreen(Unit):
    # Assume 3-deck vibrating screen
    _units = {'Area': 'ft^2'}
    _kwargs = Splitter._kwargs
    _run = Splitter._run
    
    #: Flow rate per area of screen per apareture (kg∕ft2-hr-mm)
    capacity = 0.2*907.18474
    
    #: Mesh opening (mm)
    mesh_opening = 8
    
    #: Maximum area of screen (ft^2)
    max_area = 200
    
    def _design(self):
        Area = self.ins[0].massnet/(self.capacity*self.mesh_opening)
        Design = self._results['Design'] 
        Design['N'] = N =  ceil(Area/self.max_area)
        Design['Area'] = Area/N
    