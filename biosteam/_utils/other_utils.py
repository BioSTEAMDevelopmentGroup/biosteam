# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
import biosteam as bst

__all__ = ('factor', 'checkbounds', 'approx2step', 'strtuple', 'Counter')

class Counter:
    __slots__ = ('msg', 'N', '_start')
    def __init__(self, msg=None, N=0):
        self._start = True
        self.msg = msg
        self.N = N
        
    def notify(self):
        print(f"{self.msg or 'counter'}: {self.N}")
        
    def restart(self, msg=None, N=0, notify=True):
        if notify: self.notify()
        self.msg = msg or self.msg
        self.N = N
        
    def count(self):
        self.N += 1
        
    def __repr__(self):
        return f"<Counter: msg={repr(self.msg)}, N={self.N}>"

# %% Number functions

def factor(base_units, new_units):
    if base_units == new_units: return 1
    else: return bst._Q(1, base_units).to(new_units).magnitude

def checkbounds(x, bounds):
    return bounds[0] < x < bounds[1]

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
        

