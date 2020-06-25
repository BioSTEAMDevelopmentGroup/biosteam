# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import test_utilities
from . import test_unit_operations
from . import test_sugarcane_fermentation_separations
from . import test_biorefineries

__all__ = (*test_utilities.__all__,
           *test_unit_operations.__all__,
           *test_sugarcane_fermentation_separations.__all__,
           *test_biorefineries.__all__,
           "test_biosteam")

from .test_utilities import *
from .test_unit_operations import *
from .test_biorefineries import *
from .test_sugarcane_fermentation_separations  import *

def test_biosteam():
    test_utilities()
    test_unit_operations()
    test_sugarcane_fermentation_separations()
    test_biorefineries()