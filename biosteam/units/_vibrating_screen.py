# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .decorators import cost
from .splitting import Splitter

__all__ = ('VibratingScreen',)

@cost('Area', ub=200, CE=567, cost=1010, n=0.91, BM=1.73, N='Number of screens')
class VibratingScreen(Splitter):
    # Assume 3-deck vibrating screen
    
    #: Flow rate per area of screen per apareture (kg∕ft2-hr-mm)
    capacity = 0.2*907.18474
    
    #: Mesh opening (mm)
    mesh_opening = 8
    
    _units = {'Area': 'ft^2'}
    
    def _design(self):
        self.design_results['Area'] = self.ins[0].F_mass/(self.capacity*self.mesh_opening)