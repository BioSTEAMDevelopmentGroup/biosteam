# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from .decorators import cost

__all__ = ('HammerMill',)
    
@cost('Flow rate', units='ton/hr', cost=4310, CE=541.7,
      n=0.78, ub=200., kW=6.17, BM=2.3)
class HammerMill(bst.Unit): 
    """
    Create hammer mill unit operation for reducing
    particle size.
    
    Parameters
    ----------
    ins : stream
        Solids.
    outs : stream
        Milled solids.
    
    """