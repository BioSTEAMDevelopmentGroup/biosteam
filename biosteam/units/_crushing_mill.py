# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from ._solids_separator import SolidsSeparator
from .decorators import cost

@cost('Flow rate', units='kg/hr', cost=1.5e6, CE=541.7,
      n=0.6, S=335e3, kW=2010, BM=2.3)
class CrushingMill(SolidsSeparator):
    """Create CrushingMill object.
    
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
    
