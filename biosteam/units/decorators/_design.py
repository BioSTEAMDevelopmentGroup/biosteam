# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:19:41 2019

@author: yoelr
"""

__all__ = ('design',)

def _design(self):
    D = self._results['Design']
    for i, j in self._design_basis: D[i] = j(self)

def design(name, units, func):    
    """Add design item."""
    def design_decorator(cls):
        if not cls._units: cls._units = {}
        cls._units[name] = units
        if hasattr(cls, '_design_basis'):
            cls._design_basis.append((name, func))
        else:
            cls._design_basis = [(name, func)]
            cls._design = _design
        return cls
    return design_decorator