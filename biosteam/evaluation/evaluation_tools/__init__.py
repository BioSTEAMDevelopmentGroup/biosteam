# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Mon Sep  2 04:58:34 2019

@author: yoelr
"""

from lazypkg import LazyPkg

LazyPkg(__name__, ['parameter', 'plot', 'in_parallel'])
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import parameter
from . import in_parallel

__all__ = (*parameter.__all__,
           *in_parallel.__all__)

from .parameter import *
from .in_parallel import *
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
