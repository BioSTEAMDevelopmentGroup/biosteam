# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ..decorators import cost
from thermosteam import Stream
import numpy as np
from ... import HeatUtility
from ...utils.stream_link_options import static
# from copy import copy

__all__ = ('CoolingTower',) #'CoolingTowerWithPowerDemand')

@static
@cost('Flow rate', 'Cooling water pump',
      S=557183, kW=1021, cost=283671, CE=551, n=0.8, BM=3.1)
@cost('Flow rate', 'Cooling tower',
      S=557183, kW=1598, cost=1375e3, CE=551, n=0.7, BM=1.5)
class CoolingTower(Facility):
    """Create a cooling tower that is cost based on flow rate of cooling water."""
    _units = {'Flow rate': 'kmol/hr'}
    _N_heat_utilities = 1
    _N_outs = _N_ins = 2
    evaporation = 0.01
    blowdown = 0.001
    def __init__(self, ID=''):
        thermo = HeatUtility.cooling_agents['Cooling water'].thermo
        self.makeup_water = makeup_water = Stream('cooling_tower_makeup_water',
                                                  thermo=thermo)
        loss = makeup_water.flow_proxy()
        loss.ID = 'evaporation_and_blowdown'
        super().__init__(ID, ('return_cooling_water', makeup_water),
                         ('cooling_water', loss), thermo=thermo)
        self.cooling_water_utilities = set()
        
    def _design(self):
        cwu = self.cooling_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'Cooling water':
                        cwu.add(hu)
        used = self._ins[0]
        
        #: Cooling water flow rate (kmol/hr)
        used.mol[0] = \
        self.design_results['Flow rate'] = \
        self.cooling_water = sum([i.flow for i in cwu])
        hu = self.heat_utilities[0]
        cw = hu.cooling_agents['Cooling water']
        self._outs[0].T = cw.T
        hu.ID = 'Cooling water'
        hu.cost = -self.cooling_water*cw.price_kmol
        self.makeup_water.mol[0] = self.cooling_water * (self.evaporation + self.blowdown)

CoolingTower._N_outs = CoolingTower._N_ins = 2
    
# class CoolingTowerWithPowerDemand(CoolingTower):
#     _has_power_utility = True
#     _N_heat_utilities = 1
#     cost_options = copy(CoolingTower.cost_items)
#     cost_options['Cooling tower'].kW = 0.1
#     def _cost(self):
#         super()._cost()
#         q = self._molar_flow # kmol/hr
#         hu = self.heat_utilities[0]
#         cw = hu.cooling_agents['Cooling water']
#         hu.ID = 'Cooling water'
#         hu.flow = -q
#         hu.cost = -q*cw.price_kmol