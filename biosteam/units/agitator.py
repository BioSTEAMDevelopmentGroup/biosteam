# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .auxiliary import Auxiliary
from . import design_tools as design

__all__ = ('Agitator',)

class Agitator(Auxiliary):
    __slots__ = ()
    
    def __init__(self, kW):
        super().__init__()
        self.add_power_utility(kW) # kW
        hp = kW * 1.34102
        self.baseline_purchase_costs['Agitator'] = design.compute_closed_vessel_turbine_purchase_cost(hp)