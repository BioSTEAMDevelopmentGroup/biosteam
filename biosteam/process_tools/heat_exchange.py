# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('heat_exchanger_utilities_from_units',
)
    
def heat_exchanger_utilities_from_units(units):
    """Return a list of heat utilities from all heat exchangers,
    including the condensers and boilers of distillation columns and
    flash vessel heat exchangers."""
    heat_utilities = sum([i.heat_utilities for i in units], ())
    return [i for i in heat_utilities if i.heat_exchanger]
    
