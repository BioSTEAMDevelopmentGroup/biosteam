# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:18:16 2018

@author: yoelr
"""

from .._unit import Unit
from .._exceptions import DesignError
from .decorators import cost, spec
from .metaclasses import splitter

def _ok(x, bounds):
    lb, up = bounds
    return lb < x < up
# Electricity: 16 hp / 200 ft diameter

@spec('Clarifier', 'Material', lambda M: 1.4 if M=='Steel' else 1)
@cost('Settling area', cost=2720, CE=567, exp=0.58, kW=0.00048355)
class Clarifier(Unit, metaclass=splitter):
    _units = {'Settling area': 'ft^2'}
    # Height of the clarifier tank from other designs, estimate (ft)
    height = 10
    
    # Setting the working bounds for different materials
    _bounds = {'Area steel' : (80, 8000), 'Area concrete': (8000, 125000)}
    
    def _design(self):
        # Heuristic settling area estimation
        # Settling area in ft^2 = overflow in gpm
        Design = self._results['Design']
        Design['Settling area'] = SetArea = self.outs[0].volnet *  4.4028
        # Checking to see which cost equation/material to use
        Steel_bounds, Concrete_bounds = self._bounds.values()
        if _ok(SetArea, Steel_bounds): Design['Material'] = 'Steel'
        elif _ok(SetArea, Concrete_bounds): Design['Material'] = 'Concrete'
        else: raise DesignError('Volumetric flow rate is out of working range.')
        