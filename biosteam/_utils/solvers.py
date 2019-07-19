# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
from .._exceptions import SolverError
#from .other_utils import count
import numpy as np

__all__ = ('bounded_secant', 'secant', 'wegstein', 'conditional_wegstein',
           'accelerated_secant')

def secant(f, x0, xtol, ytol, args=(), maxiter=50):
    x1 = x0*(1 + 1e-4) + (-1e-4 if x0 < 0 else 1e-4)
    pos = abs
    y0 = f(x0, *args)
    if pos(y0) < ytol: return x0
    y1 = f(x1, *args) 
    for iter in range(maxiter): 
        x1 = x0 - y1*(x1-x0)/(y1-y0)
        if (pos(x1-x0) < xtol) or (pos(y1) < ytol): return x1
        x0 = x1
        y0 = y1
        y1 = f(x1, *args)
    raise SolverError(f'failed to converge after {maxiter} iterations')

def bounded_secant(f, x0, x, x1, y0, y1, xtol, ytol, args=()):
    if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
    y = f(x, *args)
    x_ = x0
    while True:
        if abs(x-x_) < xtol: return x
        if y > ytol:
            x_ = x1 = x
            y1 = y
        elif y < ytol:
            x_ = x0 = x
            y0 = y
        else: return x
        x = x0 - y0*(x1-x0)/(y1-y0)
        y = f(x, *args)

def accelerated_secant(f, x0, xtol, ytol, args=(), maxiter=50):
    x1 = x0*(1 + 1e-4) + (-1e-4 if x0 < 0 else 1e-4)
    pos = abs
    y0 = f(x0, *args)
    if pos(y0) < ytol: return x0
    y1 = f(x1, *args)
    if pos(y1) < ytol: return x1 
    x1 = g0 = x0 - y1*(x1-x0)/(y1-y0)
    if (pos(x1-x0) < xtol): return x1
    g1 = x0
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

def wegstein(f, x0, xtol, args=(), maxiter=50):
    len_ = len(x0)
    ones = np.ones(len_)
    s = np.zeros(len_)
    x1 = y0 = f(x0, *args)
    y1 = x0
    for iter in range(maxiter):
        if (abs(y1 - x1) < xtol).all(): return x1
        y1 = f(x1, *args)
        x_diff = x1 - x0
        pos = x_diff != 0
        s[pos] = (y1[pos] - y0[pos])/x_diff[pos]
        x0 = x1
        y0 = y1
        s[s > 0.9] = 0.9
        w = ones/(ones-s)
        x1 = w*y1 + (1-w)*x1
    raise SolverError(f'failed to converge after {maxiter} iterations')

def conditional_wegstein(f, x0):
    len_ = len(x0)
    ones = np.ones(len_)
    s = np.zeros(len_)
    y0, condition = f(x0)
    x1 = y0
    while condition:
        y1, condition = f(x1)
        x_diff = x1 - x0
        pos = x_diff != 0
        s[pos] = (y1[pos] - y0[pos])/x_diff[pos]
        x0 = x1
        y0 = y1
        s[s > 0.9] = 0.9
        w = ones/(ones-s)
        x1 = w*y1 + (1-w)*x1
    return x1


    

#    if y0*y1 < 0:
#             x = x0 - y1*(x1-x0)/(y1-y0)
#             return bounded_secant(f, x0, x, x1, y0, y1, xtol, ytol, args, maxiter-iter)


