# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from flexsolve import IQ_interpolation
from typing import Callable, Optional, Tuple, Any

__all__ = ('BoundedNumericalSpecification',)

class BoundedNumericalSpecification:
    __slots__ = (
        'f', 'x0', 'x1', 'y0', 'y1', 'x', 'xtol', 'ytol', 'args', 
        'maxiter', 'checkroot', 'checkiter', 'checkbounds', 'x_last',
    )
    
    def __init__(self, 
            f: Callable,
            x0: float, 
            x1: float, 
            y0: Optional[float]=None, 
            y1: Optional[float]=None, 
            x: Optional[float]=None,
            xtol: float=0.,
            ytol: float=5e-8,
            args: Tuple[Any, ...]=(), 
            maxiter: int=50,
            checkroot: bool=False, 
            checkiter: bool=True, 
            checkbounds: bool=True
        ):
        self.f = f
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.x = x
        self.xtol = xtol
        self.ytol = ytol
        self.args = args
        self.maxiter = maxiter
        self.checkroot = checkroot
        self.checkiter = checkiter
        self.checkbounds = checkbounds
        self.x_last = None
        
    def __call__(self):
        self.x = IQ_interpolation(
            self.f, self.x0, self.x1, self.y0, self.y1, self.x, self.xtol, self.ytol, 
            self.args, self.maxiter, self.checkroot, self.checkiter, self.checkbounds,
        )
        return self.x
    
    def compile_path(self, unit): pass
    
    def create_temporary_connections(self, unit): pass