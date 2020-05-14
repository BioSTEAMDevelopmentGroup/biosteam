# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:07:01 2020

@author: yoelr
"""

__all__ = ('heat_exchanger_utilities_from_units',
)
    
def heat_exchanger_utilities_from_units(units):
    """Return a list of heat utilities from all heat exchangers,
    including the condensers and boilers of distillation columns and
    flash vessel heat exchangers."""
    heat_utilities = sum([i.heat_utilities for i in units], ())
    return [i for i in heat_utilities if i.heat_exchanger]
    
