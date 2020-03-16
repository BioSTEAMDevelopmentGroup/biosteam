# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:43:44 2020

@author: yoelr
"""
from scipy.optimize import brentq

__all__ = ('BoundedNumericalSpecification',)

class BoundedNumericalSpecification:
    __slots__ = ('f', 'a', 'b', 'kwargs')
    
    def __init__(self, f, a, b, **kwargs):
        self.f = f
        self.a = a
        self.b = b
        self.kwargs = kwargs
        
    def __call__(self):
        return brentq(self.f, self.a, self.b, **self.kwargs)
    
    def __repr__(self):
        return f"{type(self).__name__}(f={self.f}, a={self.a}, b={self.b}, kwargs={self.kwargs}"
    