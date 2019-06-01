# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from ... import Unit
from ..decorators import cost, design

__all__ = ('CoolingTower', 'CoolingTowerWithPowerDemand')

@cost('Flow rate', cost=1100, CE=567, exp=0.68)
@design('Flow rate', 'gpm', lambda self: 0.07934*self.ins[0].molnet)
class CoolingTower(Unit):
    _N_ins = 1; _N_outs = 1
    _linkedstreams = True

   
class CoolingTowerWithPowerDemand(CoolingTower):
    _has_power_utility = True
    _N_heat_utilities = 1
    power_consumption = 0.03 #: kWhr per kmol cooling water
    def _cost(self):
        super()._cost()
        q = self.ins[0].molnet # kmol/hr
        self._power_utility(q*self.power_consumption)
        hu = self._heat_utilities[0]
        cw = hu.cooling_agents['Cooling water']
        hu.duty = None
        hu.ID = 'Cooling water'
        hu.flow = -q
        hu.cost = -q*cw['Price (USD/kmol)']