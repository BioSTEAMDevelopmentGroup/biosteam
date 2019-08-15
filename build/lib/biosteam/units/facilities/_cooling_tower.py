# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ..decorators import cost
# from copy import copy

__all__ = ('CoolingTower',) #'CoolingTowerWithPowerDemand')

@cost('Flow rate', 'Cooling water pump',
      S=557183, kW=1021, cost=283671, CE=551, n=0.8, BM=3.1)
@cost('Flow rate', 'Cooling tower',
      S=557183, kW=1598, cost=1375e3, CE=551, n=0.7, BM=1.5)
class CoolingTower(Facility):
    """Create a cooling tower that is cost based on flow rate of cooling water."""
    _units = {'Flow rate': 'kmol/hr'}
    _N_heat_utilities = 1
    _N_outs = _N_ins = 0
    def __init__(self, ID=''):
        Facility.__init__(self, ID, None, None)
        self.cooling_water_utilities = set()
        
    def _design(self):
        cwu = self.cooling_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u._heat_utilities:
                    if hu.ID == 'Cooling water':
                        cwu.add(hu)
        #: Cooling water flow rate (kmol/hr)
        self._Design['Flow rate'] = self.cooling_water = q = sum([i.flow for i in cwu])
        hu = self._heat_utilities[0]
        cw = hu.cooling_agents['Cooling water']
        hu.ID = 'Cooling water'
        hu.cost = -q*cw.price_kmol
        

   
# class CoolingTowerWithPowerDemand(CoolingTower):
#     _has_power_utility = True
#     _N_heat_utilities = 1
#     cost_options = copy(CoolingTower.cost_items)
#     cost_options['Cooling tower'].kW = 0.1
#     def _cost(self):
#         super()._cost()
#         q = self._molar_flow # kmol/hr
#         hu = self._heat_utilities[0]
#         cw = hu.cooling_agents['Cooling water']
#         hu.ID = 'Cooling water'
#         hu.flow = -q
#         hu.cost = -q*cw.price_kmol