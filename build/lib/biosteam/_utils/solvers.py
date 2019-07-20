# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
from .._exceptions import SolverError
#from .other_utils import count
import numpy as np

__all__ = ('bounded_secant', 'secant', 'wegstein', 'conditional_wegstein',
           'accelerated_secant', 'accelerated_bounded_secant')

def secant(f, x0, x1, xtol, ytol=1e-7, args=(), maxiter=50):
    pos = abs
    y0 = f(x0, *args)
    if pos(y0) < ytol: return x0
    y1 = f(x1, *args) 
    if pos(y1) < ytol: return x1
    x_diff = x1-x0 
    for iter in range(maxiter): 
        x1 = x0 - y1*(x_diff)/(y1-y0)
        x_diff = x1-x0 
        if pos(x_diff) < xtol or pos(y1) < ytol: return x1
        x0 = x1
        y0 = y1
        y1 = f(x1, *args)
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_secant(f, x0, x, x1, xtol, ytol=1e-7, args=()):
    y0 = f(x0, *args)
    y1 = f(x1, *args)
    assert y0*y1 < 0, ('objective function bounderies '
                       'must have opposite signs')
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    x_ = x - 2*xtol
    while abs(x-x_) > xtol:
        y = f(x, *args)
        if y > ytol:
            x_ = x1 = x
            y1 = y
        elif y < ytol:
            x_ = x0 = x
            y0 = y
        else: return x
        x = x0 - y0*(x1-x0)/(y1-y0)
    return x

def accelerated_secant(f, x0, x1, xtol, ytol=1e-7, args=(), maxiter=50):
    pos = abs
    y0 = f(x0, *args)
    if pos(y0) < ytol: return x0
    y1 = f(x1, *args)
    if pos(y1) < ytol: return x1 
    x1 = g0 = x0 - y1*(x1-x0)/(y1-y0)
    if (pos(x1-x0) < xtol): return x1
    y0 = y1
    for iter in range(maxiter): 
        y1 = f(x1, *args)
        g1 = x1 - y1*(x1-x0)/(y1-y0)
        if (pos(g1-x1) < xtol) or (pos(y1) < ytol): return g1
        w = (x1-x0)/(x1-g1 + g0-x0)
        x1 = w*g1 + (1-w)*x1
        x0 = x1
        y0 = y1
        g0 = g1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def accelerated_bounded_secant(f, x0, x, x1, xtol, ytol=1e-7, args=(), maxiter=50):
    y0 = f(x0, *args)
    y1 = f(x1, *args)
    assert y0*y1 < 0, ('objective function bounderies '
                       'must have opposite signs')
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    
    x_old = x
    y = f(x, *args)
    if y > ytol:
        x1 = x
        y1 = y
    elif y < ytol:
        x0 = x
        y0 = y
    else: return x
    x = g0 = x0 - y0*(x1-x0)/(y1-y0)
    
    while abs(x-x_old) > xtol:
        y = f(x, *args)
        if y > ytol:
            x1 = x
            y1 = y
        elif y < ytol:
            x0 = x
            y0 = y
        else: return x
        g1 = x0 - y0*(x1-x0)/(y1-y0)
        w = (x-x_old)/(x-g1 + g0-x_old)
        x_old = x
        x = w*g1 + (1-w)*x
        if x < x0 or x > x1: x = g1
        g0 = g1
    return x

def wegstein(f, x0, xtol, args=(), maxiter=50):
    x1 = y0 = f(x0, *args)
    for iter in range(maxiter):
        if (abs(y0 - x0) < xtol).all(): return x1
        try: y1 = f(x1, *args)
        except:
            x1 = y1
            y1 = f(x1, *args)
        w = (x1-x0)/(x1-y1 + y0-x0)
        x0 = x1
        y0 = y1
        w[np.logical_not(np.isfinite(w))] = 1
        x1 = w*y1 + (1-w)*x1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def conditional_wegstein(f, x0):
    y0, condition = f(x0)
    x1 = y0
    while condition:
        try: y1, condition = f(x1)
        except:
            x1 = y1
            y1, condition = f(x1)
        w = (x1-x0)/(x1-y1 + y0-x0)
        x0 = x1
        y0 = y1
        w[np.logical_not(np.isfinite(w))] = 1
        x1 = w*y1 + (1-w)*x1
    return x1


    

#    if y0*y1 < 0:
#             x = x0 - y1*(x1-x0)/(y1-y0)
#             return bounded_secant(f, x0, x, x1, y0, y1, xtol, ytol, args, maxiter-iter)


