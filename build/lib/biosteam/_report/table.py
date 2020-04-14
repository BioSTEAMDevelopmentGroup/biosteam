# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:48:34 2018

@author: yoelr
"""
import numpy as np
import pandas as pd
from warnings import warn
import openpyxl
from .._tea import TEA, CombinedTEA
from thermosteam import Stream
from thermosteam.base import get_dimensionality, convert, stream_units_of_measure
from thermosteam.exceptions import DimensionError
import os

DataFrame = pd.DataFrame
ExcelWriter = pd.ExcelWriter

__all__ = ('stream_table', 'cost_table', 'save_system_results',
           'save_report', 'results_table', 'heat_utilities_table',
           'power_utilities_table')

def _stream_key(s):
    num = s.ID[1:]
    if num.isnumeric(): return int(num)
    else: return -1

# %% Helpful functions

def _save(tables, writer, sheet_name='Sheet1', n_row=1):
    """Save a list of tables as an excel file.
    
    Parameters
    ----------
    tables : iterable[DataFrame]
    writer : ExcelWritter
        
    """
    for t in tables:
        label = t.columns.name
        t.to_excel(writer, sheet_name, 
                   startrow=n_row, index_label=label)
        n_row += len(t.index) + 2
    
    return n_row

# %% Units

def save_report(system, file='report.xlsx', dpi='300', **stream_properties):
    """Save a system report as an xlsx file.
    
    Parameters
    ----------
    file : str
        File name to save report
    
    **stream_properties : str
        Additional stream properties and units as key-value pairs (e.g. T='degC', flow='gpm', H='kW', etc..)
        
    """
    writer = ExcelWriter(file)
    units = list(system._costunits)
    try:
        system.diagram('thorough', file='flowsheet', dpi=str(dpi), format='png')
    except:
        diagram_completed = False
        warn(RuntimeWarning('failed to generate diagram through graphviz'), stacklevel=2)
    else:
        try:
            # Assume openpyxl is used
            worksheet = writer.book.create_sheet('Flowsheet')
            flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
            worksheet.add_image(flowsheet, anchor='A1')
        except:
            # Assume xlsx writer is used
            worksheet = writer.book.add_worksheet('Flowsheet')
            flowsheet.insert_image('A1', 'flowsheet.png')
        diagram_completed = True
    
    if system.TEA:
        tea = system.TEA
        if isinstance(tea, CombinedTEA):
            costs = [cost_table(i) for i in tea.TEAs]
            _save(costs, writer, 'Itemized costs')
        else:
            # Cost table
            cost = cost_table(tea)
            cost.to_excel(writer, 'Itemized costs')
        
        # Cash flow
        tea.get_cashflow_table().to_excel(writer, 'Cash flow')
    else:
        warn(RuntimeWarning(f'Cannot find TEA object in {repr(system)}. Ignoring TEA sheets.'), stacklevel=2)
    
    
    # Stream tables
    # Organize streams by chemicals first
    streams_by_chemicals = {}
    for i in system.streams:
        if not i: continue
        chemicals = i.chemicals
        if chemicals in streams_by_chemicals:
            streams_by_chemicals[chemicals].append(i)
        else:
            streams_by_chemicals[chemicals] = [i]
    stream_tables = []
    for streams in streams_by_chemicals.values():
        stream_tables.append(stream_table(streams, **stream_properties))
    _save(stream_tables, writer, 'Stream table')
    
    # Heat utility tables
    heat_utilities = heat_utilities_table(units)
    n_row = _save(heat_utilities, writer, 'Utilities')
    
    # Power utility table
    power_utility = power_utilities_table(units)
    power_utility.to_excel(writer, 'Utilities', 
                           index_label='Electricity',
                           startrow=n_row)
    
    # General desing requirements
    results = results_table(units)
    _save(results, writer, 'Design requirements')
    writer.save()
    if diagram_completed: os.remove("flowsheet.png")

save_system_results = save_report

def results_table(units):
    """Return a list of results tables for each unit type.

    Parameters
    ----------
    units : iterable[Unit]
        
    Returns
    -------
    tables : list[DataFrame]
    
    """
    units.sort(key=(lambda u: u.line))
    
    # Organize units by units of measure:
    organized = {}
    for u in units:
        uom = (*u._units.keys(), u.line)
        if uom in organized: organized[uom].append(u)
        else: organized[uom] = [u]
    
    # Make a list of tables, keeping all results with same keys in one table
    tables = []
    for units in organized.values():
        # First table with units of measure
        table = None
        while (table is None) and units:
            u, *units = units
            table = u.results(include_utilities=False,
                              include_total_cost=False)
        if table is None: continue
        for u in units[1:]:
            table[u.ID] = u.results(with_units=False, include_utilities=False,
                                   include_total_cost=False)
        table.columns.name = (u.line, '')
        tables.append(table)
        table = u.results()
    return tables
    
def cost_table(tea):
    """Return a cost table as a pandas DataFrame object.

    Parameters
    ----------
    units : iterable[Unit]
        
    Returns
    -------
    table : DataFrame

    """
    columns = ('Unit operation',
              f'Purchase cost (10^6 USD)',
              f'Utility cost (10^6 USD/yr)')
    units = tea.units
    operating_days = tea.operating_days
    N_units = len(units)
    array = np.empty((N_units, 3), dtype=object)
    IDs = []
    types = array[0:, 0]
    C_cap = array[0:, 1]
    C_op = array[0:, 2]
    
    # Get data
    for i in range(N_units):
        unit = units[i]
        types[i] = unit.line
        C_cap[i] = unit.purchase_cost / 1e6
        C_op[i] = unit.utility_cost * operating_days * 24  / 1e6
        IDs.append(unit.ID)
    
    df = DataFrame(array, columns=columns, index=IDs)    
    if not tea.lang_factor:
        df['Installation cost (10^6 USD)'] = [u.installation_cost for u in units]
    
    return df

def heat_utilities_table(units):
    """Return a list of utility tables for each heat utility source.
    
    Parameters
    ----------
    units : iterable[Unit]
        
    Returns
    -------
    tables : list[DataFrame]
        
    """
    # Sort heat utilities by unit type, then by utility Type
    units = sorted(units, key=(lambda u: type(u).__name__))
    source = {}
    heat_utils = []
    for u in units:
        hus = u.heat_utilities
        if not hus: continue
        for i in hus: source[i] = u
        heat_utils.extend(hus)

    # Organize heatutility by ID
    heat_utils_dict = {}
    for i in heat_utils:
        ID = i.ID
        if ID in heat_utils_dict: heat_utils_dict[ID].append(i)
        else: heat_utils_dict[ID] = [i]
    
    # First table and set Type to compare with
    hu = heat_utils[0]
    Type = hu.ID
    
    # Make a list of tables, keeping all results with same Type in one table
    tables = []
    for Type, heat_utils in heat_utils_dict.items():
        data = []; index = []
        for hu in heat_utils:
            data.append((source[hu].line, hu.duty, hu.flow, hu.cost))
            index.append(source[hu].ID)
        table = DataFrame(data, index=index,
                          columns=('Unit operation',
                                   'Duty (kJ/hr)',
                                   'Flow (kmol/hr)',
                                   'Cost (USD/hr)'))
        table.columns.name = Type
        tables.append(table)
    return tables
    

def power_utilities_table(units):
    # Sort power utilities by unit type
    units = sorted(units, key=(lambda u: type(u).__name__))
    units = [u for u in units if u.power_utility]
    power_utilities = [u.power_utility for u in units]
    lenght = len(power_utilities)
    data = []
    for i, u, pu in zip(range(lenght), units, power_utilities):
        data.append((u.line, pu.rate, pu.cost))
    return DataFrame(data, index=[u.ID for u in units if u.power_utility],
                     columns=('Unit Operation', 'Rate (kW)', 'Cost (USD/hr)'))


# %% Streams

def stream_table(streams, flow='kg/hr', **props) -> 'DataFrame':
    """Return a stream table as a pandas DataFrame object.

    Parameters

    streams : array_like[Stream]
    flow : str
        Units for flow rate.
    props : str
        Additional stream properties and units as key-value pairs
    
    """
    
    # Prepare rows and columns
    ss = sorted([i for i in streams if i.ID], key=_stream_key)
    chemical_IDs = ss[0].chemicals.IDs
    n = len(ss)
    m = len(chemical_IDs)
    p = len(props)
    array = np.empty((m+p+5, n), dtype=object)
    IDs = n*[None]
    sources = array[0, :]
    sinks = array[1, :]
    phases = array[2, :]
    prop_molar_data = array[3:3+p+1,:]
    flows = array[p+3, :]
    array[p+4, :] = ''
    fracs = array[p+5:m+p+5, :]
    for j in range(n):
        s = ss[j]
        sources[j] = str(s.source or '-')
        sinks[j] = str(s.sink or '-')
        IDs[j] = s.ID
        phase = ''
        for i in s.phase:
            if i == 'l':
                phase += 'liquid|'
            elif i == 'L':
                phase += 'LIQUID|'
            elif i == 'g':
                phase += 'gas|'
            elif i == 's':
                phase += 'solid|'
        phase = phase.rstrip('|')
        phases[j] = phase
        flow_j = s.get_flow(flow)
        flows[j] = net_j = sum(flow_j)
        fracs[:,j] = flow_j/net_j if net_j > 0 else 0
        i = 0
        for attr, units in props.items():
            prop_molar_data[i, j] = s.get_property(attr, units)
            i += 1
    
    # Set the right units
    i = 0
    prop_molar_keys = [f'{attr} ({unit})' for attr, unit in props.items()]
    
    # Make data frame object
    index = ('Source', 'Sink', 'Phase')  + tuple(prop_molar_keys) + (f'flow ({flow})', 'Composition:') + tuple(chemical_IDs)
    return DataFrame(array, columns=IDs, index=index)


