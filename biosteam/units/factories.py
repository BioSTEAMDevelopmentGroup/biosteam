# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import Unit, units
from ..utils import format_unit_name
from . import decorators
from warnings import warn
import pandas as pd
import numpy as np

_add_cost = decorators.add_cost

__all__ = ('df2unit', 'xl2dct', 'xl2mod')

def df2unit(clsname, cost_items, *, supercls=None, metacls=None):
    """Return Unit subclass from cost_options DataFrame."""
    superclasses = (supercls,) if supercls else (Unit,)
    if not metacls: metacls = type(supercls)
    cls = metacls.__new__(metacls, clsname, superclasses, {})
    for ID in cost_items:
        _add_cost(cls, ID, *[(None if i is False else i) for i in cost_items[ID]], None)
    return cls

def df2dct(df):
    dct = {}
    for name_sim in df.columns.levels[0]:
        if name_sim == 'Unit': continue
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

def xl2dct(file, sheet_name=0, stacklevel=2):
    """Return dictionary of unit subclasses from excel file."""
    df_raw = pd.read_excel(file, header=[0, 1], index_col=0, nrows=12, dtype=object)
    
    # Pandas 1.2 has problems with reading excel file; columns that do not exist
    # are sometimes created with all NaN values. This block removes these columns
    # and warns when columns with real (user-specified) values are removed due 
    # to NaN values.
    df = df_raw.dropna(axis=1) 
    removed_columns = set(df_raw).difference(df)
    removed_columns = [i for i in removed_columns if not pd.isnull(df_raw[i]).all()]
    if removed_columns:
        warn(f"cost items {removed_columns} removed due to failure in reading "
              "column in excel file",
             stacklevel=stacklevel)
    index = [
        'Basis',
        'Units',
        'Size',
        'Lower bound',
        'Upper bound',
        'CEPCI',
        'Cost (USD)',
        'Exponent',
        'Electricity (kW)',
        'Installation factor',
        'Number',
        'Lifetime (yr)',
    ]
    if list(df.index) != index:
        raise RuntimeError(f"could not read '{file}'; index at column A must "
                           f"have the following entries {index}")
    return df2dct(df)

def xl2mod(file, module, sheet_name=0):
    dct = xl2dct(file, sheet_name, stacklevel=3)
    for i, j in dct.items():
        setattr(module, i, j)
        j.__module__ = module.__name__
    if not hasattr(module, '__all__'):
        module.__all__ = tuple(dct)
    
    
    
    
    
    