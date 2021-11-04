# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import Facility
from ..decorators import cost
from ... import HeatUtility
import numpy as np

__all__ = ('RefrigerationPackage',)

class RefrigerationPackage(Facility):
    """
    Create a refrigeration package that converts refrigeration usage to electricity demand [1]_.
    
    Parameters
    ----------
    ID : str, optional
        Unit ID.
    
    References
    ----------
    .. [1] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost
        Estimation". In Product and Process Design Principles: Synthesis,
        Analysis, and Evaluation (pp. 450-455). New York: Wiley.
    
    """
    ticket_name = 'RP'
    network_priority = 0
    _N_heat_utilities = 1
    _units = {'Duty': 'kJ/hr'}
    def __init__(self, ID='', agent=None, efficiency=0.35):
        raise NotImplementedError('not yet ready for users')
        self.agent = agent or HeatUtility.get_cooling_agent('propane')
        self.efficiency = efficiency
        super().__init__(ID)
        
    def _load_refrigeration_utilities(self):
        self.refrigeration_utilities = ru = set()
        ID = self.agent.ID
        for u in self.other_units:
            if u is self: continue
            for hu in u.heat_utilities:
                agent = hu.agent 
                if agent and agent.ID == ID: ru.add(hu)
        
    def _design(self):
        self._load_refrigeration_utilities()
        ru = self.refrigeration_utilities
        self.design_results['Duty'] = duty = sum([i.duty for i in ru])        
        hu, = self.heat_utilities
        hu.mix_from(ru)
        hu.reverse()
        self.power_utility.consumption = duty /3600 * self.efficiency
        