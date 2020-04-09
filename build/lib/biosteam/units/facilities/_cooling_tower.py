# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 18:28:55 2019

@author: yoelr
"""
from . import Facility
from ..decorators import cost
from ... import HeatUtility
# from copy import copy

__all__ = ('CoolingTower',) #'CoolingTowerWithPowerDemand')

@cost('Flow rate', 'Cooling water pump',
      S=557183, kW=1021, cost=283671, CE=551, n=0.8, BM=3.1)
@cost('Flow rate', 'Cooling tower',
      S=557183, kW=1598, cost=1375e3, CE=551, n=0.7, BM=1.5)
class CoolingTower(Facility):
    """
    Create a cooling tower with capital cost and power based on the flow rate 
    of cooling water as in [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    network_priority = 1
    _units = {'Flow rate': 'kmol/hr'}
    _N_heat_utilities = 1
    _N_outs = _N_ins = 2
    evaporation = 0.01
    blowdown = 0.001
    def __init__(self, ID=''):
        cooling_water = HeatUtility.get_cooling_agent('cooling_water')
        self.makeup_water = makeup_water = cooling_water.to_stream('cooling_tower_makeup_water')
        loss = makeup_water.flow_proxy()
        loss.ID = 'evaporation_and_blowdown'
        super().__init__(ID, ('return_cooling_water', makeup_water),
                         (cooling_water.to_stream(), loss), thermo=cooling_water.thermo)
        self.cooling_water_utilities = set()
        
    def _design(self):
        cwu = self.cooling_water_utilities
        if not cwu:
            for u in self.system.units:
                if u is self: continue
                for hu in u.heat_utilities:
                    if hu.ID == 'cooling_water':
                        cwu.add(hu)
        used = self._ins[0]
        
        hu = self.heat_utilities[0]
        hu.mix_from(cwu)
        
        used.imol['7732-18-5'] = \
        self.design_results['Flow rate'] = \
        self.cooling_water = hu.flow 
        
        self._outs[0].T = hu.inlet_utility_stream.T
        self.makeup_water.mol[0] = self.cooling_water * (self.evaporation + self.blowdown)
        hu.reverse()

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