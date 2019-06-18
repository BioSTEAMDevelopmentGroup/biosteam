# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 00:46:41 2019

@author: yoelr
"""
from .. import Unit, units
from . import metaclasses
from . import decorators
import pandas as pd

_finalize = decorators._extend._finalize
_index = decorators._cost._index
_cost = decorators._cost._cost
_design = decorators._design._design
design_center = decorators._design_center.design_center

__all__ = ('df2unit', 'xl2dct', 'xl2mod')

def df2unit(clsname, cost_options, *, superclass=None, metaclass=None):
    """Return Unit subclass from cost_options DataFrame."""
    dct = {'_cost': _cost,
           '_finalize': _finalize,
           'cost_options': cost_options}
    
    superclasses = (superclass,) if superclass else (Unit,)
    if not metaclass: metaclass = type(superclass)
    
    try: assert all(cost_options.index==_index), "'cost_options' index is incorrect"
    except Exception as err:
        if not isinstance(cost_options, pd.DataFrame):
            raise ValueError(f"'cost_options' must be a DataFrame, not '{type(cost_options).__name__}'")
        else:
            raise err
    
    dct['_design'] = _design
    dct['_design_basis'] = []
    add_basis = dct['_design_basis'].append
    if any(cost_options.loc['Electricity (kW)', :]):
        dct['_has_power_utility'] = True
    cls = metaclass.__new__(metaclass, clsname, superclasses, dct)
    
    done = set()
    for i in cost_options:
        column = cost_options[i]
        name = column['Basis']
        units = column['Units']
        if name in done: continue
        cls._units[name] = units
        basis = design_center(name, (units, cls._N_ins, cls._N_outs))    
        add_basis((name, basis))   
        done.add(name)
            
    return cls

def xl2dct(file, sheet_name=0):
    """Return dictionary of unit subclasses from excel file."""
    df = pd.read_excel(file, header=[0, 1])
    dct = {}
    for name_sim in df.columns.levels[0]:
        cost_options = df[name_sim]
        if '-' in name_sim:
            sim, name = name_sim.split('-')
        else:
            name = name_sim
            sim = 'static'
        name = ''.join([i.capitalize() for i in name.split(' ')])
        for i in name:
            if not (i is ' ' or i.isalnum()):
                name = name.replace(i, ' ')
        metaname = sim.casefold()
        if metaname in metaclasses.__dict__:
            metaclass = getattr(metaclasses, metaname)
            new = df2unit(name, cost_options, metaclass=metaclass)
        else:
            supername = ''.join([i.capitalize() for i in sim.split(' ')])
            if supername in units.__dict__:
                superclass = getattr(units, supername)
            else:
                raise ValueError(f"invalid simulation option '{sim}'")
            new = df2unit(name, cost_options, superclass=superclass)
        dct[name] = new
    return dct

def xl2mod(file, module, sheet_name=0):
    dct = xl2dct(file, sheet_name)
    for i, j in dct.items():
        setattr(module, i, j)
        j.__module__ = module.__name__
    if not hasattr(module, '__all__'):
        module.__all__ = tuple(dct)
    
    
    
    
    
    