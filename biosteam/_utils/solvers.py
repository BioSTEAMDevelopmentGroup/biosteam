# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
from .._exceptions import SolverError
from .other_utils import Counter
import numpy as np

__all__ = ('false_position', 'secant', 'wegstein_secant',
           'conditional_wegstein', 'aitken_secant', 'aitken',
           'wegstein', 'bounded_wegstein', 'conditional_aitken',
           'bounded_aitken', 'IQ_interpolation')

def secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    """Secant solver."""
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    dx = x1-x0 
    for iter in range(maxiter): 
        y1 = f(x1, *args)
        x1 = x0 - y1*dx/(y1-y0)
        dx = x1-x0 
        if _abs(dx) < xtol or _abs(y1) < ytol: return x1
        x0 = x1
        y0 = y1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def IQ_interpolation(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """Inverse quadratic interpolation solver."""
    _abs = abs
    if y1 < 0.: x0, y0, x1, y1 = x1, y1, x0, y0
    dx1 = dx0 = x1-x0
    f0 = yval-y0    
    if not (x0 < x < x1 or x1 < x < x0):
        try:
            # False position
            x = x0 + f0*dx0/(y1-y0)
        except ZeroDivisionError:
            # Bisection
            x = (x0+x1)/2
    yval_ub = yval + ytol
    yval_lb = yval - ytol
    while _abs(dx1) > xtol:
        y = f(x)
        if y > yval_ub:
            y2 = y1
            x2 = x1
            x1 = x
            y1 = y
        elif y < yval_lb:
            y2 = y0
            x2 = x0
            x0 = x
            y0 = y
            f0 = yval-y
        else: break
        dx1 = x1-x0
        try:
            # Inverse quadratic interpolation
            f1 = yval - y1
            f2 = yval - y2
            d01 = f0-f1
            d02 = f0-f2
            d12 = f1-f2
            f0_d12 = f0/d12
            f1_d02 = f1/d02
            f2_d01 = f2/d01
            x = x0*f1_d02*f2_d01 - x1*f0_d12*f2_d01 + x2*f0_d12*f1_d02
            assert x0 < x < x1 or x1 < x < x0
        except:
            try:
                # False position
                x = x0 + f0*dx1/(y1-y0)
            except:
                # Bisection
                x = (x0+x1)/2
            else:
                # Overshoot to prevent getting stuck
                x = x + 0.1*(x1 + x0 - 2.*x)*(dx1/dx0)**3
        dx0 = dx1
    return x

def false_position(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver."""
    _abs = abs
    if y1 < 0.: x0, y0, x1, y1 = x1, y1, x0, y0
    dx = x1-x0
    dy = yval-y0
    if not (x0 < x < x1 or x1 < x < x0):
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
        if _abs(x - x_old) < dx/10: x = (x1 + x0)/2.
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
        x0 = x1
        try:
            w = dx/(dx-g1+g0)
            x1 = w*g1 + (1.-w)*x1
        except:
            x1 = g1
        dx = x1-x0
        if _abs(dx) < xtol or _abs(y1) < ytol: return x1
        y0 = y1
        g0 = g1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_wegstein(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver with Wegstein acceleration."""
    _abs = abs
    if y1 < 0.: x0, y0, x1, y1 = x1, y1, x0, y0
    dy = yval-y0
    x_old = x = x if x0 < x < x1 or x1 < x < x0 else x0+dy*(x1-x0)/(y1-y0)
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
    else:
        return x
    dx1x0 = x1-x0
    x = g0 = x0 + dy*dx1x0/(y1-y0)
    while _abs(dx1x0) > xtol:
        y = f(x)
        if y > yval_ub:
            x1 = x
            y1 = y
        elif y < yval_lb:
            x0 = x
            y0 = y
            dy = yval - y
        else: break
        dx1x0 = x1-x0
        g1 = x0 + dy*dx1x0/(y1-y0)
        dx = x - x_old
        try:
            w = dx/(dx-g1+g0)
            x_old = x
            x = w*g1 + (1.-w)*x
        except:
            x = g0 = g1
        else:
            if x0 < x < x1 or x1 < x < x0: g0 = g1                
            else: x = g0 = g1
    return x

def wegstein(f, x0, xtol, args=(), maxiter=50):
    """Iterative Wegstein solver."""
    x1 = g0 = f(x0, *args)
    w = np.ones_like(x0)
    np_abs = np.abs
    for iter in range(maxiter):
        dx = x1-x0
        try: g1 = f(x1, *args)
        except:
            x1 = g1
            g1 = f(x1, *args)
        if (np_abs(g1-x1) < xtol).all(): return g1
        dummy = dx-g1+g0
        mask = np_abs(dummy) > 1e-16
        w[mask] = dx[mask]/dummy[mask]
        x0 = x1
        g0 = g1
        x1 = w*g1 + (1-w)*x1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def conditional_wegstein(f, x0):
    """Conditional iterative Wegstein solver."""
    g0, condition = f(x0)
    g1 = x1 = g0
    w = np.ones_like(x0)
    np_abs = np.abs
    while condition:
        try: g1, condition = f(x1)
        except:
            x1 = g1
            g1, condition = f(x1)
        dx = x1-x0
        dummy = dx-g1+g0
        mask = np_abs(dummy) > 1e-16
        w[mask] = dx[mask]/dummy[mask]
        x0 = x1
        g0 = g1
        x1 = w*g1 + (1.-w)*x1

def aitken_secant(f, x0, x1, xtol, ytol=5e-8, args=(), maxiter=50):
    """Secant solver with Aitken acceleration."""
    _abs = abs
    y0 = f(x0, *args)
    if _abs(y0) < ytol: return x0
    dx = x1-x0
    for iter in range(maxiter):
        y1 = f(x1, *args)
        if y1 == y0: return x1
        x0 = x1 - y1*dx/(y1-y0) # x0 = g
        dx = x0-x1
        if _abs(dx) < xtol or _abs(y1) < ytol: return x0
        y0 = y1
        y1 = f(x0, *args)
        if y1 == y0: return x0
        x2 = x0 - y1*dx/(y1-y0) # x2 = gg
        if _abs(dx) < xtol or _abs(y1) < ytol: return x2
        dx = x1 - x0 # x - g
        denominator = x2 + dx - x0
        x1 = x1 - dx**2./denominator if denominator else x2
        dx = x1 - x0
        y0 = y1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_aitken(f, x0, x1, y0, y1, x, yval, xtol, ytol):
    """False position solver with Aitken acceleration."""
    _abs = abs
    if y1 < 0.: x0, y0, x1, y1 = x1, y1, x0, y0
    dx1 = x1-x0
    dy = yval-y0
    if not (x0 < x < x1 or x1 < x < x0):
        x = x0 + dy*dx1/(y1-y0)
    yval_ub = yval + ytol
    yval_lb = yval - ytol
    while _abs(dx1) > xtol:
        y = f(x)
        if y > yval_ub:
            x1 = x
            y1 = y
        elif y < yval_lb:
            x0 = x
            y0 = y
            dy = yval-y
        else: 
            return x
        dx0 = x1-x0
        g = x0 + dy*dx0/(y1-y0)
        if _abs(dx0) < xtol:
            return g
        
        y = f(g)
        if y > yval_ub:
            x1 = g
            y1 = y
        elif y < yval_lb:
            x0 = g
            y0 = y
            dy = yval-y
        else:
            return g
        dx1 = x1-x0
        gg = x0 + dy*dx1/(y1-y0)
        dxg = x - g
        try: x = x - dxg**2./(gg + dxg - g)
        except:
            # Add overshoot to prevent getting stuck
            x = gg + 0.1*(x1+x0-2*gg)*(dx1/dx0)**3. 
        else:
            if not (x0 < x < x1 or x1 < x < x0):
                x = gg + 0.1*(x1+x0-2*gg)*(dx1/dx0)**3. 
    return x

def aitken(f, x, xtol, args=(), maxiter=50):
    """Iterative Aitken solver."""
    gg = x
    x = x.copy()
    abs_ = abs
    np_abs = np.abs
    for iter in range(maxiter):
        try: g = f(x, *args)
        except:
            x = gg.copy()
            g = f(x, *args)
        dxg = x - g
        if (abs_(dxg) < xtol).all(): return g
        gg = f(g, *args)
        dgg_g = gg - g
        if (abs_(dgg_g) < xtol).all(): return gg
        dummy = dgg_g + dxg
        mask = np_abs(dummy) > 1e-16
        x[mask] -= dxg[mask]**2/dummy[mask]
    raise SolverError(f'failed to converge after {maxiter} iterations')
    
def conditional_aitken(f, x):
    """Conditional iterative Aitken solver."""
    condition = True
    x = x.copy()
    gg = x
    np_abs = np.abs
    while condition:
        try:
            g, condition = f(x)
        except:
            x = gg.copy()
            g, condition = f(x)
        if not condition: return g
        gg, condition = f(g)
        dxg = x - g
        dummy = gg + dxg - g
        mask = np_abs(dummy) > 1e-16
        x[mask] -= dxg[mask]**2/dummy[mask]