# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
cost = bst.decorators.cost

__all__ = ('CIPpackage',)


# %%
@cost('Flow rate', units='kg/hr',
      S=63, cost=421e3, CE=522, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    ticket_name = 'CIP'
    line = 'CIP Package'
    network_priority = 0
    _N_ins = 1
    _N_outs = 1
    