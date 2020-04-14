# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:18:16 2018

@author: yoelr
"""

from ._splitter import Splitter
from ..exceptions import DesignError
from .decorators import cost

__all__ = ('Clarifier',)

_iswithin = lambda x, bounds: bounds[0] < x < bounds[1]
# Electricity: 16 hp / 200 ft diameter

@cost('Settling area', cost=2720, CE=567, n=0.58, kW=0.0005)
class Clarifier(Splitter):
    _units = {'Settling area': 'ft^2'}
    # Height of the clarifier tank from other designs, estimate (ft)
    height = 10
    
    # Setting the working bounds for different materials
    _bounds = {'Area steel' : (80, 8000), 'Area concrete': (8000, 125000)}
    
    def _design(self):
        # Heuristic settling area estimation
        # Settling area in ft^2 = overflow in gpm
        Design = self.design_results
        Design['Settling area'] = SetArea = self.outs[0].F_vol *  4.4028
        # Checking to see which cost equation/material to use
        Steel_bounds, Concrete_bounds = self._bounds.values()
        if _iswithin(SetArea, Steel_bounds): Design['Material'] = 'Steel'
        elif _iswithin(SetArea, Concrete_bounds): Design['Material'] = 'Concrete'
        else: raise DesignError('Volumetric flow rate is out of working range.')
        

_cost = Clarifier._cost
def _extended_cost(self):
    _cost(self)
    self.purchase_costs['Clarifier'] *= 1.4 if self.design_results['Material']=='Steel' else 1
Clarifier._cost = _extended_cost
        
        
        