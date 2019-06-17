# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:05:53 2019

@author: yoelr
"""
import biosteam as bst
import pandas as pd
from copy import copy
from ._extend import extend_finalize
from math import ceil

__all__ = ('cost',)

_index = pd.Index(['Basis',
                   'Units',
                   'Size',
                   'Limit',
                   'CEPCI',
                   'Cost (USD)',
                   'Exponent',
                   'Electricity (kW)'])

def _cost(self):
    table = self.cost_options
    Cost = self._results['Cost']
    Design = self._results['Design']
    net_e = 0
    for i in table:
        basis, units, value, limit, CE, cost, n, e = table[i]
        size = Design[basis]
        if limit:
            Design['#'+i] = N = ceil(size/limit)
            f = size/value
            S = f/N
            Cost[i] = N*bst.CEPCI/CE*cost*S**n
            net_e += e*f
        else:
            S = size/value
            Cost[i] = bst.CEPCI/CE*cost*S**n
            net_e += e*S
    if net_e: self._power_utility(net_e)

def cost(basis, name=None, *, cost, exp, CE, S=1, kW=0, limit=None):    
    r"""Add item purchase cost based on exponential scale up:
    
    :math:`C_{f.o.b.} = (N)(CE)(cost)(\frac{basis}{S})^{exp}` 
    
    :math:`Electricity\ rate = (N)(kW)(\frac{basis}{S})`
    
    Where `N` corresponds to the minimum number of parallel units given the size limit (if any).
    
    **Parameters**
    
        **basis:** Name of Design item used for scaling.
        
        **name:** Name of item.
        
        **cost:** Cost of item.
        
        **exp:** Exponential factor.
        
        **CE:** Chemical engineering plant cost index.
        
        **S:** Size.
        
        **kW:** Electricity rate.
        
        **limit:** Size limit.
        
    **Examples**
    
        :doc:`Unit decorators`
    
    """
    def cost_decorator(cls):
        extend_finalize(cls)
        if kW: cls._has_power_utility = True
        try:
            units = cls._units[basis]
        except KeyError:
            raise RuntimeError(f'units of cost basis ({basis}) is not available in "{cls.__name__}._units" dictionary')
        data = [basis, units, S, limit, CE, cost, exp, kW]
        if hasattr(cls, 'cost_options'):
            if 'cost_options' not in cls.__dict__:
                cls.cost_options = copy(cls.cost_options)
            if not name:
                raise ValueError("must pass a 'name' for purchase cost item")
            if name in cls.cost_options:
                raise ValueError("name '{name_}' not available, must pass a name not previously used")
            cls.cost_options[name or cls.line] = data
        else:
            if '_cost' in cls.__dict__:
                raise RuntimeError("cannot decorate class as '_cost' method is already implemented")
            columns = pd.Index((name or cls.line,), name='Item')
            cls.cost_options = pd.DataFrame(data, _index, columns)
            cls._cost = _cost
        return cls
    return cost_decorator

