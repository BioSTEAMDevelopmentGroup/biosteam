# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Mon Mar  4 16:40:50 2019

@author: yoelr
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
"""
from .decorators import cost
from ..utils.unit_warnings import lb_warning
from .._unit import Unit

__all__ = ('ScrewFeeder',)

@cost('Flow rate', ub=10e4, CE=567, cost=1096, n=0.22)
class ScrewFeeder(Unit):
    length = 30 #: ft
    _N_outs = 1
    _has_power_utility = True
    _minimum_flow = 400
    _units = {'Flow rate': 'ft^3/hr'}
    
    def _design(self):
        feed = self.ins[0]
        r = self.results
        Design = r['Design']
        F_vol = feed.F_vol*35.315 # ft3/hr
        if F_vol < self._minimum_flow:
            lb_warning('Flow rate', F_vol, self._minimum_flow)
        Design['Flow rate'] = F_vol
<<<<<<< HEAD
        F_mass = feed.F_mass*0.0006124 #lb/s
=======
        F_mass = 0.0006124 * feed.F_mass # lb/s
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
        self._power_utility(0.0146*F_mass**0.85*self.length*0.7457)