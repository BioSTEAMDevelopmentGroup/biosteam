# -*- coding: utf-8 -*-
<<<<<<< HEAD
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
from thermosteam import exceptions
from thermosteam.exceptions import *

__all__ = ('DesignError', *exceptions.__all__)
del exceptions

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""
<<<<<<< HEAD

=======
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
