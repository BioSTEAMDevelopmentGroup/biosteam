# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Thu Aug 23 22:14:01 2018

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
from ._solids_separator import SolidsSeparator
from .decorators import cost

__all__ = ('CrushingMill',)

@cost('Flow rate', units='kg/hr', cost=1.5e6, CE=541.7,
      n=0.6, S=335e3, kW=2010, BM=2.3)
class CrushingMill(SolidsSeparator):
<<<<<<< HEAD
    """Create CrushingMill object.
=======
    """
    Create CrushingMill object.
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    
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
    
