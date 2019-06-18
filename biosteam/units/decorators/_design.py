# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:19:41 2019

@author: yoelr
"""
from ._design_center import design_center
__all__ = ('design',)

def _design(self):
    D = self._results['Design']
    for i, j in self._design_basis: D[i] = j(self)

def design(name, units, func=None):    
    """Add design item.
    
    **Parameters**
    
        **name:** Name of design item.
        
        **units:** Units of measure of design item.
        
        **func:** Should return design item given the Unit object. If None, defaults to function predefined for given name and units.
        
    **Examples**
    
        :doc:`Unit decorators`
    
    """
    def design_decorator(cls):
        f = func or design_center(name, (units, cls._N_ins, cls._N_outs))
        if not cls._units:
            cls._units = {}
        elif '_units' not in cls.__dict__:
            cls._units = cls._units.copy()
        cls._units[name] = units
        if hasattr(cls, '_design_basis'):
            cls._design_basis.append((name, f))
        else:
            cls._design_basis = [(name, f)]
            cls._design = _design
        return cls
    return design_decorator

