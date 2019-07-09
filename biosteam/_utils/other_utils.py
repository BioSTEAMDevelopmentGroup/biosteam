# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
from biosteam import _Q

__all__ = ('factor', 'checkbounds', 'approx2step', 'strtuple', 'isbetween', 'count')

def count():
    i = 0
    print('starting count')
    while True:
        i += 1
        print(i)
        yield

def isbetween(x0, x, x1):
    if x0 > x1: x0, x1 = x1, x0
    return x0 < x < x1


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

def strtuple(iterable):
    """Return string of all items in the tuple""" 
    string = ''
    function = type(strtuple)
    for i in iterable:
        if isinstance(i , function):
            string += i.__name__ + ', '
        else:
            string += str(i) + ', '
    string = string.rstrip(', ')
    string = '(' + string + ')'
    return string
        

