# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biorefineries.tests.test_biorefineries import (
    test_oilcane_O6, 
    test_oilcane_O8, 
    test_oilcane_O9,
    test_sugarcane,
    test_cornstover,
)
    
if __name__ == '__main__':
    test_oilcane_O6()
    test_oilcane_O8() 
    test_oilcane_O9()
    test_sugarcane()
    test_cornstover()