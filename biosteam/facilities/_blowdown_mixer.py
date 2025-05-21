# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam._graphics import mixer_graphics
import biosteam as bst

__all__ = ('BlowdownMixer',)

class BlowdownMixer(bst.Facility, bst.Mixer):
    network_priority = 2
    _graphics = mixer_graphics
    ticket_name = 'BDM'