# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Mon Mar  4 11:32:01 2019

@author: yoelr
"""

=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
from .._unit import Unit
from .decorators import cost

__all__ = ('MagneticSeparator',)

@cost('Flow rate', units='kg/hr', CE=576, cost=533471, S=333333, n=0.6)
class MagneticSeparator(Unit): pass
    