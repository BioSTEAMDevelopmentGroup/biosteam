# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from ... import Unit
from ..decorators import cost, design
from ..metaclasses import static

__all__ = ('CoolingTower', 'CoolingTowerWithPowerDemand')

@cost('Flow rate', cost=1100, CE=567, exp=0.7)
@design('Flow rate', 'gpm')
class CoolingTower(Unit, metaclass=static):
    """Create a cooling tower that is cost based on flow rate of cooling water.
    
    **ins**
    
        [0] Cooling water
        
    **outs**
    
        [0] Cooling water
    
    """

   
class CoolingTowerWithPowerDemand(CoolingTower):
    _has_power_utility = True
    _N_heat_utilities = 1
    power_consumption = 0.020 #: kWhr per kmol cooling water
    def _cost(self):
        super()._cost()
        q = self.ins[0].molnet # kmol/hr
        self._power_utility(q*self.power_consumption)
        hu = self._heat_utilities[0]
        cw = hu.cooling_agents['Cooling water']
        hu.ID = 'Cooling water'
        hu.flow = -q
        hu.cost = -q*cw['Price (USD/kmol)']