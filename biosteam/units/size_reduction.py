# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for size reduction.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.Shredder

.. autoclass:: biosteam.units.HammerMill

"""
from .._unit import Unit
from .decorators import cost

__all__ = ('Shredder', 'HammerMill')

@cost('Flow rate', units='kg/hr', cost=2.5e6,
      CE=567.3, n=0.6, S=500e3, kW=3000, BM=1.39)
class Shredder(Unit):  pass

# Default kW based on an industrial corn dry-grind hammer mill
@cost('Flow rate', units='ton/hr', cost=4310, lb=2, ub=200,
      CE=567, n=0.78, kW=6.17, BM=2.3) 
class HammerMill(Unit):  pass