# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
from biosteam import _Q

__all__ = ('factor', 'checkbounds', 'approx2step', 'strtuple', 'function',
           'golden_ratio')

# %% Solvers

def golden_ratio(f, x0, xf, xtol, guess):
    e0 = f(x0)
    ef = f(xf)
    guesses = (guess-xtol, guess+xtol)
    for x in guesses:
        e = f(x)
        if e > 0:
            xf = x
            ef = e
        else:
            x0 = x
            e0 = e
    if x0 in guesses and xf in guesses: return guess
    elif ef < 0: raise ValueError('f(xf) must be positive')
    elif e0 > 0: raise ValueError('f(x0) must be negative')
    gr_v = 0.618
    gr_h = 0.382
    if abs(ef) > abs(e0):
        x = gr_v*xf + gr_h*x0
    else:
        x = gr_v*x0 + gr_h*xf
    dx = xtol+1
    while dx > xtol:
        e = f(x)
        if e > 0:
            dx = xf - x
            xf = x
            ef = e
        else:
            dx = x - x0
            x0 = x
            e0 = e
        if ef > -e0:
            x = gr_v*xf + gr_h*x0
        else:
            x = gr_v*x0 + gr_h*xf
    return x


# %% Number functions

def factor(base_units, new_units):
    if base_units == new_units: return 1
    else: return _Q(1, base_units).to(new_units).magnitude

def checkbounds(x, bounds):
    lb, up = bounds
    return lb < x < up

def approx2step(val, x0, dx):
    """Approximate value, val, to closest increment/step, dx, starting from x0."""
    while True:
        if x0 > val: break
        x0 += dx
    return x0


# %% String functions

function = type(checkbounds)

def strtuple(iterable):
    """Return string of all items in the tuple""" 
    string = ''
    for i in iterable:
        if isinstance(i , function):
            string += i.__name__ + ', '
        else:
            string += str(i) + ', '
    string = string.rstrip(', ')
    string = '(' + string + ')'
    return string
        

