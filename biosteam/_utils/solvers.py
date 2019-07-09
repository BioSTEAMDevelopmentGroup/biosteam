# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
from .._exceptions import SolverError
import numpy as np

__all__ = ('boundsolve', 'wegstein')

def boundsolve(f, x0, x, x1, xtol, ytol, args=(), maxiter=50):
    y0 = f(x0, *args)
    y1 = f(x1, *args)
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    y = f(x)
    x_ = x0
    it = 0
    while abs(x-x_) > xtol:
        it += 1
        if it > maxiter:
            raise SolverError('failed to converge')
        if y > ytol:
            x_ = x1 = x
            y1 = y
        elif y < ytol:
            x_ = x0 = x
            y0 = y
        else:
            break
        dx = x1-x0
        x = x0 + - y0*(dx)/(y1-y0)
        y = f(x, *args)
    return x, y

def wegstein(f, x, xtol, args=(), maxiter=50):
    # Prepare variables
    len_ = len(x)
    ones = np.ones(len_)
    s = np.zeros(len_)

    # First run
    x0 = x
    x1 = gx0 = f(x, *args)
    gx1 = x0

    # Check convergence
    abs_ = abs

    # Outer loop
    it = 0
    while (abs_(gx1 - x1) > xtol).all():
        it += 1
        if it > maxiter: raise SolverError('failed to converge')
        # Get relaxation factor and set next iteration
        gx1 = f(x1, *args)
        x_diff = x1 - x0
        pos = x_diff != 0
        s[pos] = (gx1[pos] - gx0[pos])/x_diff[pos]
        x0 = x1
        gx0 = gx1
        s[s > 0.9] = 0.9
        w = ones/(ones-s)
        x1 = w*gx1 + (1-w)*x1
    return x1