# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
import doctest

__all__ = ('test_heat_utility',
           'test_power_utility',
           'test_utilities',
)

def test_heat_utility():
    from biosteam import _heat_utility
    doctest.testmod(_heat_utility)
    
def test_power_utility():
    from biosteam import _power_utility
    doctest.testmod(_power_utility)
    
def test_utilities():
    test_heat_utility()
    test_power_utility()
    
if __name__ == '__main__':
    test_utilities()