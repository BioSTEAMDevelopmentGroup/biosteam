# -*- coding: utf-8 -*-
"""
Created on Wed May 22 23:14:02 2019

@author: yoelr
"""
from .. import _ureg
from .._exceptions import DimensionError

__all__ = ('DisplayUnits',)

setattr_ = object.__setattr__
getattr_ = object.__getattribute__

class DisplayUnits:
    
    def __init__(self, **display_units):
        setattr_(self, '_display_units', display_units)
        list_keys = []
        for k, v in display_units.items():
            try:
                dims = getattr(_ureg, v).dimensionality
            except:
                try: 
                    dims = [getattr(_ureg, i).dimensionality for i in v]
                    list_keys.append(k)
                except: 
                    dims = v
            setattr_(self, k, dims)
        for k in list_keys:
            display_units[k] = display_units[k][0]
    
    def __setattr__(self, name, unit):
        if name not in self.__dict__:
            raise AttributeError(f"can't set display units for '{name}'")
        if not isinstance(unit, str) and isinstance(unit, type(self._display_units[name])):
            self._display_units[name] = unit
            return
        name_dim = getattr_(self, name)
        unit_dim = getattr(_ureg, unit).dimensionality
        if isinstance(name_dim, list):
            if unit_dim not in name_dim:
                raise DimensionError(f"dimensions for '{name}' must be either {', '.join(name_dim[:-1])} or {name_dim[-1]}; not {unit_dim}")    
        else:
            if name_dim != unit_dim:
                raise DimensionError(f"dimensions for '{name}' must be in {name_dim}, not {unit_dim}")
        self._display_units[name] = unit
    
    def __getattribute__(self, name):
        if name[0] == '_':
            return getattr_(self, name)
        else: 
            return self._display_units[name]
    
    def __iter__(self):
        for i in self._display_units.values():
            if isinstance(i, list): yield i[0]
            yield i
            
    def __repr__(self):
        sig = ', '.join((f"{i}='{j}'" if isinstance(j, str) else f'{i}={j}') for i,j in self._display_units.items())
        return f'{type(self).__name__}({sig})'
        