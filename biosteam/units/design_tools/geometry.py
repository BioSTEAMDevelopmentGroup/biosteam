# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from flexsolve import njitable
from math import pi

__all__ = ('cylinder_diameter_from_volume', 
           'cylinder_volume_from_diameter')

@njitable
def cylinder_diameter_from_volume(volume, length_to_diameter):
    return (4. * length_to_diameter * volume / pi)**(1./3.)

@njitable
def cylinder_volume_from_diameter(diameter, length_to_diameter):
    return pi * (diameter / 2) ** 2 * diameter / length_to_diameter