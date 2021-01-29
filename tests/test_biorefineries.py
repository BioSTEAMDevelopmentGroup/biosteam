# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biorefineries.tests.test_biorefineries import *


def ignore_import_error(f):
    def g():
        try: f()  
        except ImportError: pass
    return g

# To run theses tests, you'll need to get the github version.
# If you have the PyPI installation, these tests will not run.
# This is OK, the loss in coverage is very small.
try: test_LAOs = ignore_import_error(test_LAOs)
except NameError: pass
try: test_lactic = ignore_import_error(test_lactic)
except NameError: pass
try: test_HP_cellulosic = ignore_import_error(test_HP_cellulosic)
except NameError: pass
try: test_HP_sugarcane = ignore_import_error(test_HP_sugarcane)
except NameError: pass
try: test_ethanol_adipic = ignore_import_error(test_ethanol_adipic)
except NameError: pass
try: test_animal_bedding = ignore_import_error(test_animal_bedding)
except NameError: pass
try: test_wheatstraw = ignore_import_error(test_wheatstraw)
except NameError: pass