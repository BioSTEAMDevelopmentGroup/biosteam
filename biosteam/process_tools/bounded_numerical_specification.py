# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from flexsolve import IQ_interpolation

__all__ = ('BoundedNumericalSpecification',)

class BoundedNumericalSpecification: # pragma: no cover
    __slots__ = ('args', 'solver', 'kwargs')
    
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        
    def __call__(self):
        return IQ_interpolation(*self.args, **self.kwargs)
    