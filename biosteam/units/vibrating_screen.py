# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 20:51:22 2019

@author: yoelr
"""
from ..unit import Unit
from .splitter import Splitter
from numpy import ceil

class VibratingScreen(Unit):
    # Assume 3-deck vibrating screen
    kwargs = Splitter.kwargs
    _run = Splitter._run
    
    #: Flow rate per area of screen per apareture (kg∕ft2-hr-mm)
    capacity = 0.2*907.18474
    
    #: Mesh opening (mm)
    mesh_opening = 8
    
    #: Maximum area of screen (ft^2)
    max_area = 200
    
    def _cost(self):
        results = self.results
        Design = results['Design']
        
        Cost = results['Cost']
        Area = self.ins[0].massnet/(self.capacity*self.mesh_opening)
        Design['N_screens'] = N_screens =  ceil(Area/self.max_area)
        Area /= N_screens
        Design['Area'] = Area
        Cost['Vibrating screens'] = N_screens * self.CEPCI/567 * 1010 * Area **0.91
    