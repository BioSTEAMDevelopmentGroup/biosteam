# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:18:16 2018

@author: yoelr
"""

from biosteam.unit import Unit
from biosteam.exceptions import DesignError
from . import Splitter
import numpy as np

def checkbounds(x, bounds):
    lb, up = bounds
    return lb < x < up

class Clarifier(Unit):
    _has_power_utility = True
    _kwargs = Splitter._kwargs
    _run = Splitter._run
    
    # Height of the clarifier tank from other designs, estimate (ft)
    height = 10
    
    # Setting the working bounds for different materials
    bounds = {'Area steel' : (80, 8000), 'Area concrete': (8000, 125000)}
    _results_UnitsOfMeasure = {'Settling area': 'ft^2',
                               'Cost': 'USD'}
    
    
    def _cost(self):
        results = self.results
        Design = results['Design']
        overflow = self.ins[0]
        
        # Heuristic settling area estimation
        # Settling area in ft^2 = overflow in gpm
        SetArea = overflow.volnet *  4.4028
        Design['Settling area'] = SetArea
        power = self._calc_energy(SetArea)*0.7457 # in kW
        self.power_utility(power) 
        
        # Checking to see which cost equation/material to use
        Steel_bounds, Concrete_bounds = self.bounds.values()
        if checkbounds(SetArea, Steel_bounds):
            Cost = 3810*SetArea**0.58
            Design['Material'] = 'Steel'
        elif checkbounds(SetArea, Concrete_bounds):
            Cost = 2720*SetArea**0.58
            Design['Material'] = 'Concrete'
        else:
            raise DesignError('Volumetric flow rate is out of working range')
        results['Cost']['Clarifier'] = Cost*self.CEPCI/567
        
    def _calc_energy(self, SetArea):
        # Finding diameter of the tank in ft
        diameter = np.sqrt(SetArea*4/np.pi)
        # Energy assuming quadratic relationship in hp
        Energy = 16 * (diameter/200)**2
        return Energy