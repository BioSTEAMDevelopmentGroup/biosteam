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
from ._splitter import Splitter

__all__ = ('VibratingScreen',)

@cost('Area', units='ft^2', ub=200, CE=567, cost=1010, n=0.91, BM=1.73, N='Number of screens',
      fsize=lambda self, _: self.ins[0].F_mass/(self.capacity*self.mesh_opening))
class VibratingScreen(Splitter):
    # Assume 3-deck vibrating screen
    
    #: Flow rate per area of screen per apareture (kg∕ft2-hr-mm)
    capacity = 0.2*907.18474
    
    #: Mesh opening (mm)
    mesh_opening = 8