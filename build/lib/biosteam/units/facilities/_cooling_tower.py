# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from ... import Unit
from ..decorators import cost
from ..metaclasses import static
from copy import copy

__all__ = ('CoolingTower', 'CoolingTowerWithPowerDemand')

@cost('Flow rate', 'Cooling tower', cost=1100, CE=567, n=0.7)
class CoolingTower(Unit, metaclass=static):
    """Create a cooling tower that is cost based on flow rate of cooling water.
    
    **ins**
    
        [0] Cooling water
        
    **outs**
    
        [0] Cooling water
    
    """
    _units = {'Flow rate': 'gpm'}
    def _design(self):
        self._Design['Flow rate'] = self._ins[0].molnet*0.07932

   
class CoolingTowerWithPowerDemand(CoolingTower):
    _has_power_utility = True
    _N_heat_utilities = 1
    cost_options = copy(CoolingTower.cost_items)
    cost_options['Cooling tower'].kW = 0.25215
    def _cost(self):
        super()._cost()
        q = self.ins[0].molnet # kmol/hr
        hu = self._heat_utilities[0]
        cw = hu.cooling_agents['Cooling water']
        hu.ID = 'Cooling water'
        hu.flow = -q
        hu.cost = -q*cw.price_kmol