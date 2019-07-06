# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
from biosteam import _Q

__all__ = ('factor', 'checkbounds', 'approx2step', 'strtuple', 'function',
           'ThermoSolver')

def count():
    i = 0
    while True:
        i += 1
        print(i)
        yield

# %% Solvers
        
class ThermoSolver:
    __slots__ = ('T', 'P', 'Q', 'V')
    
    tolerance = {'T': 0.001,
                 'P': 100,
                 'Q': 0.2,
                 'V': 0.001}
    
    def __init__(self):
        self.T = self.P = self.Q = self.V = None
    
    def __call__(self, xvar, yvar, f, x0, x1, y0, y1, yval):
        xtol = self.tolerance[xvar]
        ytol = self.tolerance[yvar]
        x = getattr(self, xvar)
        if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
        if x is None or not (x0 < x < x1 ):
            x = x0 + (yval-y0)*(x1-x0)/(y1-y0)
        dx = xtol+1
        y = ytol+1
        while abs(dx) > xtol and abs(yval-y) > ytol:
            y = f(x)
            if y > yval:
                dx = x1 - x
                x1 = x
                y1 = y
            else:
                dx = x - x0
                x0 = x
                y0 = y
            x = x0 + (yval-y0)*(x1-x0)/(y1-y0)
        setattr(self, xvar, x)
        setattr(self, yvar, y)
        return x

    def __repr__(self):
        return f"<{type(self).__name__}: T={self.T}, P={self.P}, Q={self.T}, V={self.V}>"


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
        

