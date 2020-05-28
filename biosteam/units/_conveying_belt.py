# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .decorators import cost
from .._unit import Unit

__all__ = ('ConveyingBelt',)

@cost('Flow rate', CE=567, cost=813, ub=2000, n=0.38, N='Number of conveyors')
class ConveyingBelt(Unit):
    length = 40 #: ft
    height = 20 #: ft
    _N_outs = 1
    _has_power_utility = True
    _minimum_flow = 120
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        self.design_results['Flow rate'] = F_vol = feed.F_vol*35.315 # ft3/hr
        if F_vol < self._minimum_flow:
            self._lb_warning('Flow rate', F_vol, self._minimum_flow)
        F_mass = feed.F_mass*0.0006124 #lb/s
        self.power_utility(0.00058*F_mass**0.82*self.length
                           + self.height*0.00182*F_mass * 0.7457) # kW
        
        

