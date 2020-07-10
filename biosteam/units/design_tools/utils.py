# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 13:54:54 2020

@author: yoelr
"""
from flexsolve import njitable

__all__ = ('approx2step',)

<<<<<<< HEAD
@njitable
=======
@njitable(cache=True)
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
def approx2step(val, x0, dx):
    """Approximate value, val, to closest increment/step, dx, starting from x0."""
    while True:
        if x0 > val: break
        x0 += dx
<<<<<<< HEAD
    return x0
=======
    return x0
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
