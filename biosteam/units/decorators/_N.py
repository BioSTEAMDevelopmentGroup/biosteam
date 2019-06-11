# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 21:06:29 2019

@author: yoelr
"""
from ... import Unit
from copy import copy
from math import ceil
from pandas import DataFrame
from ._extend import extend_summary
__all__ =('N',)

_index = ('Basis', 'Limit')

def _N(self):
    table = self.N_options
    Design = self._results['Design']
    for i in table:
        basis, ub = table[i]
        value = Design[basis]
        Design[i] = N = ceil(value/ub)
        Design[basis] = value/N

def N(basis, name='N', *, limit):
    """Add design item for number of parallel unit items based on basis limit.
    
    **Parameters**
    
        **basis:** Name of Design item used for scaling.
        
        **name:** Name of item.
        
        **limit:** Maximum value of basis.
        
    **Examples**
    
        :doc:`Unit decorators`
    
    """
    def N_decorator(cls):
        extend_summary(cls)
        if '_N' in cls.__dict__:
            raise RuntimeError("cannot decorate class as '_N' method is already implemented")
        if hasattr(cls, 'N_options'):
            if 'N_options' not in cls.__dict__:
                cls.N_options = copy(cls.N_options)
            if name in cls.N_options:
                raise ValueError("name '{name_}' not available, must pass a name not previously used")
            cls.N_options[name] = (basis, limit)
        else:
            cls.N_options = DataFrame([basis, limit],
                                      columns=(name,),
                                      index=_index)
            cls._N = _N
        return cls
    return N_decorator