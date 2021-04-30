# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 13:54:54 2020

@author: yoelr
"""
from numba import njit

__all__ = ('approx2step',)

@njit(cache=True)
def approx2step(val, x0, dx):
    """Approximate value, val, to closest increment/step, dx, starting from x0."""
    while True:
        if x0 > val: break
        x0 += dx
    return x0
