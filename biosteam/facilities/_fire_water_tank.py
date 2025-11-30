# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
cost = bst.decorators.cost

__all__ = ('FireWaterTank',)

@cost('Flow rate', 'Tank', S=8343, units='kg/hr',
      CE=522, cost=803000, n=0.7, BM=1.8)
@cost('Flow rate', 'Pump', S=8343, units='kg/hr',
      CE=522, cost=15000, n=0.8, BM=1.7, kW=94.3375)
class FireWaterTank(bst.Facility):
    ticket_name = 'FWT'
    network_priority = 0
    _N_ins = _N_outs = 0
    def _init(self, fire_water_flow_rate=8343):
        self.fire_water_flow_rate = fire_water_flow_rate # Assumes 4 hr storage capacity
    
    def _design(self):
        self.design_results['Flow rate'] = self.fire_water_flow_rate