# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .._unit import Unit
from .decorators import cost

__all__ = ('Shredder',)

@cost('Flow rate', units='kg/hr', cost=2.5e6,
      CE=567.3, n=0.6, S=500e3, kW=3000, BM=1.39)
class Shredder(Unit):  pass