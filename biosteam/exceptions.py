# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import exceptions
from thermosteam.exceptions import *

__all__ = ('DesignError', *exceptions.__all__)
del exceptions

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""
    
class GraphicsWarning(RuntimeWarning):
    """RuntimeWarning regarding diagrams."""
    
class FailedEvaluation(RuntimeWarning):
    """RuntimeWarning regarding failed model evaluation."""
