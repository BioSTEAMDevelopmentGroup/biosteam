# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
def test_sugarcane_biorefinery():
    from biorefineries import sugarcane as sc
    from biorefineries.tests.test_biorefineries import test_sugarcane
    sc._system_loaded = False
    test_sugarcane()

def test_cornstover_biorefinery():
    from biorefineries import cornstover as cs
    from biorefineries.tests.test_biorefineries import test_cornstover
    cs.Biorefinery.cache.clear()
    test_cornstover()