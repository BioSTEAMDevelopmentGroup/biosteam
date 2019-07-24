# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
from .._exceptions import SolverError
# from .other_utils import count
import numpy as np

__all__ = ('bounded_secant', 'secant', 'wegstein_secant',
           'conditional_wegstein', 'aitken_secant', 'aitken',
           'wegstein', 'bounded_wegstein', 'conditional_aitken')

def secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args) 
    if _abs(y1) < ytol: return x1
    x_diff = x1-x0 
    for iter in range(maxiter): 
        x1 = x0 - y1*(x_diff)/(y1-y0)
        x_diff = x1-x0 
        if _abs(x_diff) < xtol or _abs(y1) < ytol: return x1
        x0 = x1
        y0 = y1
        y1 = f(x1, *args)
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_secant(f, x0, x, x1, xtol, ytol=1e-7, args=()):
    y0 = f(x0, *args)
    _abs = abs
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args)
    if _abs(y1) < ytol: return x1
    assert y0*y1 < 0, ('objective function bounderies '
                       'must have opposite signs')
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    delx = x1-x0
    x = x0 - y0*delx/(y1-y0)
    while _abs(delx) > xtol:
        y = f(x, *args)
        if y > ytol:
            x1 = x
            y1 = y
        elif y < -ytol:
            x0 = x
            y0 = y  
        else: break
        delx = x1-x0
        x = x0 - y0*delx/(y1-y0)
    return x

def wegstein_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args)
    if _abs(y1) < ytol: return x0
    g0 = x1 - y1*(x1-x0)/(y1-y0)
    y0 = y1
    delx = g0-x1
    x1 = g0
    for iter in range(maxiter):
        y1 = f(x1, *args)
        g1 = x1 - y1*delx/(y1-y0)
        w = delx/(delx-g1+g0)
        x0 = x1
        x1 = w*g1 + (1-w)*x1
        delx = x1-x0
        if (_abs(delx) < xtol) or (_abs(y1) < ytol): return x1
        y0 = y1
        g0 = g1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_wegstein(f, x0, x, x1, xtol, ytol=1e-7, args=(), maxiter=50):
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args)
    if _abs(y1) < ytol: return x1
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
    delx = x-x_old
    delx_bounds = x1-x0
    while _abs(delx_bounds) > xtol:
        y = f(x, *args)
        if y > ytol:
            x1 = x
            y1 = y
        elif y < -ytol:
            x0 = x
            y0 = y
        else: return x
        delx_bounds = x1-x0
        g1 = x0 - y0*delx_bounds/(y1-y0)
        w = (delx)/(delx-g1+g0)
        x_old = x
        x = w*g1 + (1-w)*x
        g0 = g1
    return x

def wegstein(f, x0, xtol, args=(), maxiter=50):
    x1 = g0 = f(x0, *args)
    logical_not = np.logical_not
    isfinite = np.isfinite
    for iter in range(maxiter):
        delx = x1-x0
        try: g1 = f(x1, *args)
        except:
            x1 = g1
            g1 = f(x1, *args)
        if (abs(g1-x1) < xtol).all(): return g1
        w = delx/(delx-g1+g0)
        w[logical_not(isfinite(w))] = 1
        x0 = x1
        g0 = g1
        x1 = w*g1 + (1-w)*x1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def conditional_wegstein(f, x0):
    g0, condition = f(x0)
    x1 = g0
    logical_not = np.logical_not
    isfinite = np.isfinite
    while condition:
        try: g1, condition = f(x1)
        except:
            x1 = g1
            g1, condition = f(x1)
        delx = x1-x0
        w = delx/(delx-g1+g0)
        x0 = x1
        g0 = g1
        w[logical_not(isfinite(w))] = 1
        x1 = w*g1 + (1-w)*x1
    return g0

def aitken_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args)
    if _abs(y1) < ytol: return x1 
    g = x1 - y1*(x1-x0)/(y1-y0)
    y0 = y1
    y1 = f(g, *args)
    if _abs(y1) < ytol: return g
    x = gg = g - y1*(g-x1)/(y1-y0)
    for iter in range(maxiter): 
        try:
            y1 = f(x, *args)
        except:
            x = gg
            y1 = f(x, *args)
        gg = x - y1*(x-g)/(y1-y0)
        if (_abs(gg-g) < xtol) or (_abs(y1) < ytol): return gg
        x = x0 - (g - x0)**2/(gg - 2*g + x0)
        x0 = x1
        y0 = y1
        g = gg
    raise SolverError(f'failed to converge after {maxiter} iterations')

def aitken(f, x, xtol, args=(), maxiter=50):
    gg = x
    abs_ = abs
    for iter in range(maxiter):
        try: g = f(x, *args)
        except:
            x = gg
            g = f(x, *args)
        if (abs_(g-x) < xtol).all(): return gg
        gg = f(g, *args)
        x = x - (g - x)**2/(gg-2*g+x)
        if (abs_(gg-g) < xtol).all(): return gg
    raise SolverError(f'failed to converge after {maxiter} iterations')
    
def conditional_aitken(f, x):
    logical_not = np.logical_not
    isfinite = np.isfinite
    condition = True
    gg = x
    while condition:
        try:
            g, condition = f(x)
        except:
            x = gg
            g, condition = f(x)
        if not condition: g
        gg, condition = f(g)
        x = x - (g - x)**2/(gg-2*g+x)
        pos = logical_not(isfinite(x))
        x[pos] = gg[pos]
    return x



# def aitken_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
#     _abs = abs
#     y0 = f(x0, *args)
#     if _abs(y0) < ytol: return x0
#     gg = x0
#     for iter in range(maxiter):
#         import pdb
#         pdb.set_trace()
#         y1 = f(x1, *args)
#         g = x1 - y1*(x1-gg)/(y1-y0)
#         if _abs(gg-g) < xtol or _abs(y1) < ytol: return g
#         y0 = f(g, *args)
#         gg = g - y0*(g-x1)/(y0-y1)
#         if _abs(gg-g) < xtol or _abs(y0) < ytol: return gg
#         x1 = x1 - (g - x1)**2/(gg-2*g+x1)
#     raise SolverError(f'failed to converge after {maxiter} iterations')