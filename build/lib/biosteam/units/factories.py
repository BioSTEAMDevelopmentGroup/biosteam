# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 00:46:41 2019

@author: yoelr
"""
from .. import Unit, units
from ..utils import format_unit_name, static_flow_and_phase
from . import decorators
import pandas as pd
import numpy as np

_add_cost = decorators.add_cost
_index = np.array(('Basis',
                   'Units',
                   'Size',
                   'Upper bound',
                   'CEPCI',
                   'Cost (USD)',
                   'Exponent',
                   'Electricity (kW)',
                   'Installation factor',
                   'Number'))


__all__ = ('df2unit', 'xl2dct', 'xl2mod')

def df2unit(clsname, cost_items, *, supercls=None, metacls=None):
    """Return Unit subclass from cost_options DataFrame."""
    superclasses = (supercls,) if supercls else (Unit,)
    if not metacls: metacls = type(supercls)
    cls = metacls.__new__(metacls, clsname, superclasses, {})
    for ID in cost_items:
        _add_cost(cls, ID, *cost_items[ID], None)
    return cls

def df2dct(df):
    dct = {}
    for name_sim in df.columns.levels[0]:
        cost_items = df[name_sim]
        if '-' in name_sim:
            sim, name = name_sim.split('-')
            is_static = False
        else:
            name = name_sim
            sim = 'Unit'
            is_static = True
        name = format_unit_name(name)
        try:
            if is_static:
                supercls = Unit
            else:
                supercls = getattr(units, sim)
        except:
            supername = ''.join([i.capitalize() for i in sim.split(' ')])
            try: supercls = getattr(units, supername)
            except AttributeError:
                raise ValueError(f"invalid simulation option '{sim}'")
        dct[name] = df2unit(name, cost_items, supercls=supercls)
    return dct

def xl2dct(file, sheet_name=0):
    """Return dictionary of unit subclasses from excel file."""
    return df2dct(pd.read_excel(file, header=[0, 1]))

def xl2mod(file, module, sheet_name=0):
    dct = xl2dct(file, sheet_name)
    for i, j in dct.items():
        setattr(module, i, j)
        j.__module__ = module.__name__
    if not hasattr(module, '__all__'):
        module.__all__ = tuple(dct)
    
    
    
    
    
    