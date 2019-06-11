# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:05:53 2019

@author: yoelr
"""
from inspect import signature
from functools import reduce
from ._extend import extend_finalize

__all__ = ('spec',)

# %% Tracking cost factors
_mul = lambda x, y: x*y

def _spec(self):
    D = self._results['Design']
    C = self._results['Cost']
    S = self._specdata
    for i in S: C[i] *= reduce(_mul, (f(D[p]) for f, p in S[i]))

def spec(item, param, func):
    if not isinstance(param, str):
        raise ValueError(f"param must be a string, not a '{type(param).__name__}' object.")
    N = len(signature(func).parameters)
    if N != 1:
        raise ValueError(f"one and only one argument in 'func' signature is allowed ({N or 'no'} parameters given)")
    
    def spec_decorator(cls):
        extend_finalize(cls)
        if hasattr(cls, '_specdata'): data = cls._specdata
        else: cls._specdata = data = {}
        if item in data: data[item].append((func, param))
        else: data[item] = [(func, param)]
        cls._spec = _spec
        return cls
    return spec_decorator
        
    
    