# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('NotImplementedMethod',)

class NotImplementedMethodType:
    __slots__ = ()
    
    @property
    def __name__(self): return "NotImplementedMethod"
    @property
    def __doc__(self): return None
    def __new__(self): return NotImplementedMethod
    def __call__(self): return NotImplemented
    def __bool__(self): return False
    def __repr__(self): return "NotImplementedMethod"

NotImplementedMethod = object.__new__(NotImplementedMethodType)