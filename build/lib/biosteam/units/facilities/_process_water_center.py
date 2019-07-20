# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 18:43:39 2019

@author: yoelr
"""
from . import Facility
import biosteam as bst

__all__ = ('ProcessWaterCenter',)

cost = bst.units.decorators.cost

@cost('Makeup water flow rate', 'Makeup water pump',
      CE=551, kW=20*0.7457, cost=6864, S=155564, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Process water pump',
      CE=551, kW=75*0.7457, cost=15292, S=518924, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Tank',
      CE=522, cost=250e3, S=451555, n=0.7, BM=1.7)
class ProcessWaterCenter(Facility):
    _N_ins = 2
    _N_outs = 1
    makeup_water_price = 0.0002113 #: (USD/kg)
    boiler_blowdown = 0.002
    RO_rejection = 0.25
    CT_evaporation = 0.001
    CT_blowdown = 0.01
    _units = {'Makeup water flow rate': 'kg/hr',
              'Process water flow rate': 'kg/hr'}
    def __init__(self, ID, ins=None, outs=(), *, BT, CT):
        Facility.__init__(self, ID, ins, outs)
        makeup_water = bst.Stream('Makeup_water', species=bst.Species('Water',)) 
        if self._ins[1]:
            self._ins[1].link = makeup_water
        else:
            self._ins[1] = makeup_water
        self.BT = BT
        self.CT = CT
        
    def _run(self):
        process_water = self._outs[0].molnet
        recycle_water = self._ins[0].molnet
        boiler = self.BT.total_steam * self.boiler_blowdown * 1/(1-self.RO_rejection)
        cooling_tower = self.CT.cooling_water * (self.CT_evaporation + self.CT_blowdown)
        makeup_water_stream = self._ins[1] 
        Design = self._Design
        process_loss = process_water - recycle_water
        if process_loss > 0:
            makeup_water = boiler + cooling_tower + process_loss
        else:
            makeup_water = boiler + cooling_tower
        makeup_water_stream.mol[0] = makeup_water
        
        Design['Process water flow rate'] = (process_water + boiler + cooling_tower)*18.015
        Design['Makeup water flow rate'] = makeup_water*18.015
        makeup_water_stream.price = self.makeup_water_price
        
        