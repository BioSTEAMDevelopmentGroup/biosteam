# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
"""
__all__ = ('NotImplementedMethod',)

class NotImplementedMethodType:
    __slots__ = ()
    
    @property
    def __name__(self): return "NotImplementedMethod"
    @property
    def __doc__(self): return None
    def __new__(self): return NotImplementedMethod
<<<<<<< HEAD
    def __call__(self): pass
=======
    def __call__(self): return NotImplemented
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
    def __bool__(self): return False
    def __repr__(self): return "NotImplementedMethod"

NotImplementedMethod = object.__new__(NotImplementedMethodType)