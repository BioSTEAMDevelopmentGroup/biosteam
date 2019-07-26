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
           'wegstein', 'bounded_wegstein', 'conditional_aitken',
           'bounded_aitken', 'isbetween')

def isbetween(x0, x, x1):
    return x0 < x < x1 if x1 > x0 else x1 < x < x0

def secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    """Secant solver."""
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

def bounded_secant(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver."""
    _abs = abs
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    dx = x1-x0
    dy = yval-y0
    if not isbetween(x0, x, x1):
        x = x0 + dy*dx/(y1-y0)
    yval_ub = yval + ytol
    yval_lb = yval - ytol
    x_old = x
    while _abs(dx) > xtol:
        y = f(x)
        if y > yval_ub:
            x1 = x
            y1 = y
        elif y < yval_lb:
            x0 = x
            y0 = y
            dy = yval-y
        else: break
        dx = x1-x0
        x = x0 + dy*dx/(y1-y0)
        if (x - x_old) < dx/3: x = (x1 + x0)/2
    return x

def wegstein_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    """Secant solver with Wegstein acceleration."""
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    y1 = f(x1, *args)
    if _abs(y1) < ytol: return x0
    g0 = x1 - y1*(x1-x0)/(y1-y0)
    y0 = y1
    dx = g0-x1
    x1 = g0
    for iter in range(maxiter):
        y1 = f(x1, *args)
        g1 = x1 - y1*dx/(y1-y0)
        w = dx/(dx-g1+g0)
        x0 = x1
        x1 = w*g1 + (1-w)*x1
        dx = x1-x0
        if (_abs(dx) < xtol) or (_abs(y1) < ytol): return x1
        y0 = y1
        g0 = g1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_wegstein(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver with Wegstein acceleration."""
    _abs = abs
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    dy = yval-y0
    x_old = x = x if isbetween(x0, x, x1) else x0+dy*(x1-x0)/(y1-y0)
    y = f(x)
    yval_ub = yval + ytol
    yval_lb = yval - ytol
    if y > yval_ub:
        x1 = x
        y1 = y
    elif y < yval_lb:
        x0 = x
        y0 = y
        dy = yval - y
    else: return x
    dx_bounds = x1-x0
    x = g0 = x0 + dy*dx_bounds/(y1-y0)
    while _abs(dx_bounds) > xtol:
        y = f(x)
        if y > yval_ub:
            x1 = x
            y1 = y
        elif y < yval_lb:
            x0 = x
            y0 = y
            dy = yval - y
        else: return x
        dx_bounds = x1-x0
        g1 = x0 + dy*dx_bounds/(y1-y0)
        dx = x - x_old
        w = (dx)/(dx-g1+g0)
        x_old = x
        x = w*g1 + (1-w)*x
        g0 = g1
        if not isbetween(x0, x, x1): x = g1
    return x

def wegstein(f, x0, xtol, args=(), maxiter=50):
    """Iterative Wegstein solver."""
    x1 = g0 = f(x0, *args)
    logical_not = np.logical_not
    isfinite = np.isfinite
    for iter in range(maxiter):
        dx = x1-x0
        try: g1 = f(x1, *args)
        except:
            x1 = g1
            g1 = f(x1, *args)
        if (abs(g1-x1) < xtol).all(): return g1
        w = dx/(dx-g1+g0)
        w[logical_not(isfinite(w))] = 1
        x0 = x1
        g0 = g1
        x1 = w*g1 + (1-w)*x1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def conditional_wegstein(f, x0):
    """Conditional iterative Wegstein solver."""
    g0, condition = f(x0)
    g1 = x1 = g0
    logical_not = np.logical_not
    isfinite = np.isfinite
    while condition:
        try: g1, condition = f(x1)
        except:
            x1 = g1
            g1, condition = f(x1)
        dx = x1-x0
        w = dx/(dx-g1+g0)
        x0 = x1
        g0 = g1
        w[logical_not(isfinite(w))] = 1
        x1 = w*g1 + (1-w)*x1

def aitken_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    """Secant solver with Aitken acceleration."""
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    dx = x1-x0
    for iter in range(maxiter):
        y1 = f(x1, *args)
        x0 = x1 - y1*dx/(y1-y0) # x0 = g
        dx = x0-x1
        if (_abs(dx) < xtol) or (_abs(y1) < ytol): return x0
        y0 = y1
        y1 = f(x0, *args)
        x2 = x0 - y1*dx/(y1-y0) # x2 = gg
        if (_abs(dx) < xtol) or (_abs(y1) < ytol): return x2
        dx = x1 - x0 # x - g
        x1 = x1 - dx**2/(x2 + dx - x0)
        dx = x1 - x0
        y0 = y1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_aitken(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver with Aitken acceleration."""
    _abs = abs
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    dx_bounds = x1-x0
    dy = yval-y0
    if not isbetween(x0, x, x1):
        x = x0 + dy*dx_bounds/(y1-y0)
    yval_ub = yval + ytol
    yval_lb = yval - ytol
    while _abs(dx_bounds) > xtol:
        y = f(x)
        if y > yval_ub:
            x1 = x
            y1 = y
        elif y < yval_lb:
            x0 = x
            y0 = y
            dy = yval-y
        else: return x
        dx_bounds = x1-x0
        g = x0 + dy*dx_bounds/(y1-y0)
        if _abs(dx_bounds) < xtol: return g
        
        y = f(g)
        if y > yval_ub:
            x1 = g
            y1 = y
        elif y < yval_lb:
            x0 = g
            y0 = y
            dy = yval-y
        else: return g
        dx_bounds = x1-x0
        gg = x0 + dy*dx_bounds/(y1-y0)
        dxg = x - g
        x = x - dxg**2/(gg + dxg - g)
        if not isbetween(x0, x, x1): x = gg
    return x

def aitken(f, x, xtol, args=(), maxiter=50):
    """Iterative Aitken solver."""
    gg = x
    abs_ = abs
    logical_not = np.logical_not
    isfinite = np.isfinite
    for iter in range(maxiter):
        try: g = f(x, *args)
        except:
            x = gg
            g = f(x, *args)
        dxg = x - g
        if (abs_(dxg) < xtol).all(): return g
        gg = f(g, *args)
        dgg_g = gg - g
        if (abs_(dgg_g) < xtol).all(): return gg
        x = x - dxg**2/(dgg_g + dxg)
        pos = logical_not(isfinite(x))
        x[pos] = gg[pos]
    raise SolverError(f'failed to converge after {maxiter} iterations')
    
def conditional_aitken(f, x):
    """Conditional iterative Aitken solver."""
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
        if not condition: return g
        gg, condition = f(g)
        dxg = x - g
        x = x - dxg**2/(gg + dxg - g)
        pos = logical_not(isfinite(x))
        x[pos] = gg[pos]