# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

__all__ = ('ethanol_subsystem',)

from ._ethanol_subsystem_example import create_ethanol_subsystem_example


def __getattr__(name):
    if name == 'ethanol_subsystem':
        global ethanol_subsystem
        ethanol_subsystem = create_ethanol_subsystem_example()
        return ethanol_subsystem
    else:
        raise AttributeError("module %s has no attribute %s" %(__name__, name))