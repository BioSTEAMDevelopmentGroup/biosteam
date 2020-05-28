# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import biorefinery_tests
from . import unit_operation_tests
from . import utility_tests

__all__ = (*biorefinery_tests.__all__,
           "test_biosteam")

from .biorefinery_tests import *
from .unit_operation_tests import *
from .utility_tests import *

def test_biosteam():
    utility_tests.test_utilities()
    unit_operation_tests.test_unit_operations()
    biorefinery_tests.test_biorefineries()