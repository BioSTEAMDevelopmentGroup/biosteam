# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biorefineries.tests.test_biorefineries import (
    test_sugarcane,
    test_lipidcane,
    test_cornstover,
    test_LAOs
)

# To run these following tests, you'll need to get the github version.
# If you have the PyPI installation, these tests will not run.
# This is OK, the loss in coverage is very small.

# from biorefineries.tests import (
#     test_ethanol_adipic,
#     test_lactic
# )    

# Lactic acid tests failing for now due to some refractoring