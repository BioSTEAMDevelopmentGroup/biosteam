# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from scipy.optimize import brentq

__all__ = ('BoundedNumericalSpecification',)

class BoundedNumericalSpecification: # pragma: no cover
    __slots__ = ('args', 'solver', 'kwargs')
    
    def __init__(self, *args, solver=brentq, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.solver = solver
        
    def __call__(self):
        return self.solver(*self.args, **self.kwargs)
    