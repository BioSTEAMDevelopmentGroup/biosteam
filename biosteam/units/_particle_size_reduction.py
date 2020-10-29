# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from ._solids_separator import SolidsSeparator
from .decorators import cost

__all__ = ('CrushingMill', 'HammerMill')

@cost('Flow rate', units='kg/hr', cost=1.5e6, CE=541.7,
      n=0.6, S=335e3, kW=2010, BM=2.3)
class CrushingMill(SolidsSeparator):
    """
    Create crushing mill unit operation for the 
    separation of sugarcane juice from the baggasse.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] Shredded sugar cane
        * [1] Recycle water
    outs : stream sequence 
        * [0] Bagasse
        * [1] Juice
    moisture_content : float
                       Fraction of water in Baggasse.
    
    """
    
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