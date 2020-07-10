# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 22:18:16 2018

@author: yoelr
"""

=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
from ._splitter import Splitter
from ..exceptions import DesignError
from .decorators import cost

__all__ = ('Clarifier',)

<<<<<<< HEAD
_iswithin = lambda x, bounds: bounds[0] < x < bounds[1]
=======
_iswithin = lambda x, bounds: bounds[0] <= x <= bounds[1]
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
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
<<<<<<< HEAD
        if _iswithin(SetArea, Steel_bounds): Design['Material'] = 'Steel'
=======
        if SetArea <= Steel_bounds[-1]: Design['Material'] = 'Steel'
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
        elif _iswithin(SetArea, Concrete_bounds): Design['Material'] = 'Concrete'
        else: raise DesignError('Volumetric flow rate is out of working range.')
        

_cost = Clarifier._cost
def _extended_cost(self):
    _cost(self)
    self.purchase_costs['Clarifier'] *= 1.4 if self.design_results['Material']=='Steel' else 1
Clarifier._cost = _extended_cost
        
        
        