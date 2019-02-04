# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:15:20 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam.units import Splitter
import numpy as np

class RotaryVacuumFilter(Unit):
    #: Revolutions per min
    rpm = 20/60
    
    #: For crystals (lb/day-ft^2)
    filter_rate = 6000
    
    kwargs = {'split': None}
    bounds = {'Individual area': (10, 800)}
    
    @staticmethod
    def _calc_Area(flow:'kg/hr', filter_rate:'lb/day-ft^2') -> 'ft^2':
        flow *= 52.91
        return flow/filter_rate
        
    _run = Splitter._run
    
    def _design(self):
        """
        'Area': (ft^2)
        'Individual area': (ft^2)
        """
        flow = sum(stream.massnet for stream in self.outs)
        area = self._calc_Area(flow, self.filter_rate)
        self.results['Design']['Area'] = area
        
    def _cost(self):
        """
        'Cost': (USD)
        """
        results = self.results
        Design = results['Design']
        Area = Design['Area']
        ub = results.bounds['Individual area'][1]
        N_vessels = np.ceil(Area/ub)
        iArea = Area/N_vessels # individual vessel
        Design['Individual area'] = iArea
        
        Cost = np.exp(11.796-0.1905*np.log(iArea)+0.0554*(np.log(iArea))**2)
        results['Cost']['Cost'] = N_vessels*Cost
        
        
RVF = RotaryVacuumFilter


