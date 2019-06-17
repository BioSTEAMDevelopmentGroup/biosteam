# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 00:46:41 2019

@author: yoelr
"""
from .. import _Q, Unit, units
from .._stream import mol_flow_dim, mass_flow_dim, vol_flow_dim
from .._exceptions import DimensionError
from . import metaclasses
from . import decorators
from types import ModuleType
import pandas as pd
import ntpath
import re

_finalize = decorators._extend._finalize
_index = decorators._cost._index
_cost = decorators._cost._cost
_design = decorators._design._design

__all__ = ('from_excel',)

flow_rates = {'kmol/hr': lambda self: self._ins[0].molnet,
              'kg/hr': lambda self: self._ins[0].massnet,
              'm3/hr': lambda self: self._ins[0].volnet}

duties = {'kJ/hr': lambda self: self._outs[0].H - self._ins[0].H}

def get_flow(units):
    q = _Q(1, units)
    dim = q.dimensionality
    if dim == mol_flow_dim:
        factor = 1/q.to('kmol/hr').magnitude
        return lambda self: self._ins[0].molnet * factor
    elif dim == mass_flow_dim:
        factor = 1/q.to('kg/hr').magnitude
        return lambda self: self._ins[0].massnet * factor
    elif dim == vol_flow_dim:
        factor = 1/q.to('m3/hr').magnitude
        return lambda self: self._ins[0].volnet * factor
    else:
        raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")

def get_duty(units):
    factor = _Q(1, 'kJ/hr').to(units)
    return lambda self: factor*(self._outs[0].H - self._ins[0].H) 

def from_options(clsname, cost_options, *, superclass=None, metaclass=None):
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
            print(cost_options)
            print(_index)
            raise err
    
    dct['_design'] = _design
    dct['_design_basis'] = []
    add_basis = dct['_design_basis'].append
    done = set()
    for i in cost_options:
        column = cost_options[i]
        basis = column['Basis']
        units = column['Units']
        basis_units = basis, units
        if basis_units in done: continue
        done.add(basis_units)
        if basis == 'Flow rate':
            if units in flow_rates:
                add_basis((basis, flow_rates[units]))
            else:
                add_basis((basis, get_flow(units)))
        elif basis == 'Duty':
            if units in duties:
                add_basis((basis, duties[units]))
            else:
                add_basis((basis, get_duty(units)))
        else:
            raise ValueError(f"unknown basis '{basis}', basis must be either 'Flow rate' or 'Duty'")
          
    return metaclass.__new__(metaclass, clsname, superclasses, dct)

def from_excel(file, module_name=None):
    df = pd.read_excel(file, header=[0, 1, 2])
    modname = module_name or re.sub(
            r"\B([A-Z])", r"_\1", ntpath.basename(file).split('.')[-2])
    module = ModuleType(modname)
    if modname in units.__dict__:
        raise ValueError(f"module name '{modname}' already in biosteam.units")
    for name_sim in zip(*df.columns.levels[:-1]):
        cost_options = df[name_sim]
        name, sim = name_sim
        name = ''.join([i.capitalize() for i in name.split(' ')])
        metaname = sim.casefold()
        if metaname in metaclasses.__dict__:
            metaclass = getattr(metaclasses, metaname)
            new = from_options(name, cost_options, metaclass=metaclass)
        else:
            supername = ''.join([i.capitalize() for i in sim.split(' ')])
            if supername in units.__dict__:
                superclass = getattr(units, supername)
            else:
                raise ValueError(f"invalid simulation option '{sim}'")
            new = from_options(name, cost_options, superclass=superclass)
        setattr(module, name, new)
    setattr(units, modname, module)
    return module
