# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from .decorators import cost

__all__ = ('MagneticSeparator',)

@cost('Flow rate', units='kg/hr', CE=576, cost=533471, S=333333, n=0.6)
class MagneticSeparator(Unit): pass
    