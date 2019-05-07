# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:05:53 2019

@author: yoelr
"""
import pandas as pd

__all__ = ('cost',)

_index = pd.MultiIndex.from_tuples(
        (('Basis', 'param'),
         ('Basis', 'units'),
         ('Basis', 'value'),
         ('Number', ''),
         ('CEPCI', ''),
         ('Cost (USD)', ''),
         ('Exponent', ''),
         ('Electricity (kW)', '')))

def _cost(self):
    table = self.cost_options
    Cost = self._results['Cost']
    Design = self._results['Design']
    check = isinstance
    net_e = 0
    for i in table:
        param, units, value, N, CEPCI, cost, n, e = table[i]
        S = Design[param]/value
        if check(N, str): N = Design[N]
        Cost[i] = N*self.CEPCI/CEPCI*cost*S**n
        net_e += (N*e*S)
    if net_e: self._power_utility(net_e)

def cost(basis, name=None, *, cost, exp, CE, S=1, kW=0, N=1):    
    """Add item purchase cost based on exponential scale up."""
    def cost_decorator(cls):
        if kW: cls._has_power_utility = True
        try: units = cls._units[basis]
        except KeyError: raise RuntimeError(f'Units of cost basis ({basis}) is not available in "{cls.__name__}._units" dictionary.')
        data = [basis, units, S, N, CE, cost, exp, kW]
        if hasattr(cls, 'cost_options'):
            if not name: raise ValueError("Must pass a 'name' for purchase cost item.")
            if name in cls.cost_options: raise ValueError("name '{name_}' in not available, must pass a name not previously used.")
            cls.cost_options[name or cls.line] = data
        else:
            if '_cost' in cls.__dict__: raise RuntimeError("Cannot decorate class as '_cost' method is already implemented.")
            columns = pd.Index((name or cls.line,), name='Item')
            cls.cost_options = pd.DataFrame(data, _index, columns)
            cls._cost = _cost
        return cls
    return cost_decorator

