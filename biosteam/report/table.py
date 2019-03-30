# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:48:34 2018

@author: yoelr
"""
from ..stream import Stream, mol_flow_dim, mass_flow_dim, vol_flow_dim
from .. import Q_, pd, np
from ..exceptions import DimensionError

DataFrame = pd.DataFrame
ExcelWriter = pd.ExcelWriter

__all__ = ('stream_table', 'cost_table', 'save_system_results',
           'results_table', 'heat_utilities_table', 'power_utilities_table')

# %% Helpful functions

def _units_wt_cost(units):
    """Remove units that do not have capital or operating cost."""
    r = 0
    units = list(units)
    for i in range(len(units)):
        i -= r
        unit = units[i]
        Summary = unit.results['Summary']
        ci = Summary['Purchase cost']
        co = Summary['Utility cost']
        if not (ci or co):
            del units[i]
            r += 1
    return units

def _units_sort_by_cost(units):
    """Sort by first grouping units by line and then order groups by maximum capital."""
    unit_lines = list({u.line for u in units})
    line_dict = {ut:[] for ut in unit_lines}
    for u in units:
        line_dict[u.line].append(u)
    units = []
    unit_lines.sort(key=lambda ut: max(u.results['Summary']['Purchase cost'] for u in line_dict[ut]), reverse=True)
    for key in unit_lines:
        ulist = line_dict[key]
        ulist.sort(key=lambda x: x.results['Summary']['Purchase cost'], reverse=True)
        units += ulist
    return units

def save(tables, writer, sheet_name='Sheet1', n_row=1):
    """Save a list of tables as an excel file.
    
    **Parameters**
    
        **tables:** iterable[DataFrame]
        
        **writer:** [ExcelWritter]
        
    """
    for t in tables:
        label = t.columns.name
        t.to_excel(writer, sheet_name, 
                   startrow=n_row, index_label=label)
        n_row += len(t.index) + 2
    
    return n_row

# %% Units

def save_system_results(system, file='report.xlsx', **stream_properties):
    """Save a system table as an xlsx file."""
    writer = ExcelWriter(file)
    units = system.units
    
    # Cost table
    cost = cost_table(system)
    cost.to_excel(writer, 'Itemized costs', 
                  index_label=cost.columns.name)
    
    # Stream tables
    table = stream_table(system.streams, **stream_properties)
    table.to_excel(writer, 'Stream table', 
                    index_label=cost.columns.name)
    
    # Heat utility tables
    heat_utilities = heat_utilities_table(units)
    n_row = save(heat_utilities, writer, 'Utilities')
    
    # Power utility table
    power_utility = power_utilities_table(units)
    power_utility.to_excel(writer, 'Utilities', 
                           index_label='Electricity',
                           startrow=n_row)
    
    # General desing requirements
    results = results_table(units)
    save(results, writer, 'Design requirements')
    writer.save()
    

def results_table(units):
    """Return a list of results tables for each unit type.

    **Parameters**

        **units:** iterable[Unit]
        
    **Returns:**
    
        **tables:** list[DataFrame]
    
    """
    units = _units_wt_cost(units)
    units.sort(key=(lambda u: u.line))
    
    # First table and set keys to compare with
    r = units[0]._results
    keys = tuple(r.nested_keys())
    table = units[0].results()
    del units[0]
    
    # Make a list of tables, keeping all results with same keys in one table
    tables = []
    for u in units:
        r = u._results
        new_keys = tuple(r.nested_keys())
        if new_keys == keys:
            table = pd.concat((table, u.results(with_units=False)), axis=1)
        else:
            tables.append(table)
            table = r.table()
            keys = new_keys
    
    # Add the last table
    tables.append(table)
    
    return tables
    
def cost_table(system):
    """Return a cost table as a pandas DataFrame object.

    **Parameters**:

        **units:** array_like[Unit]
         
    **Returns**
    
        **unit_table:** [DataFrame] Units as indexes with the following columns
            * 'Unit Type': Type of unit
            * 'CI (10^6 USD)': Capital investment
            * 'UC (10^6 USD/yr)': Annual utility cost

    """
    columns = ('Type',
              f'Fixed Capital Investment (10^6 USD)',
              f'Utility Cost (10^6 USD/yr)')
    units = system.units
    units = _units_wt_cost(units)
    units = _units_sort_by_cost(units)
    
    # Initialize data
    try:
        o = system.tea.options
    except AttributeError:
        if not hasattr(system, 'tea'):
            raise ValueError('Cannot find TEA object related to system. A TEA object of the system must be created first.')
    lang_factor = o['Lang factor']
    operating_days = o['Operating days']
    N_units = len(units)
    array = np.empty((N_units, 3), dtype=object)
    IDs = []
    types = array[0:, 0]
    C_cap = array[0:, 1]
    C_op = array[0:, 2]
    
    # Get data
    for i in range(N_units):
        unit = units[i]
        Summary = unit.results['Summary']
        types[i] = unit.line
        C_cap[i] = Summary['Purchase cost'] * lang_factor / 1e6
        C_op[i] = Summary['Utility cost'] * operating_days * 24  / 1e6
        IDs.append(unit.ID)
    
    return DataFrame(array, columns=columns, index=IDs)    

def heat_utilities_table(units):
    """Return a list of utility tables for each heat utility source.
    
    **Parameters**
    
        **units:** iterable[units]
        
    **Returns**
    
        **tables:** list[DataFrame]
        
    """
    # Sort heat utilities by unit type, then by utility Type
    units = sorted(units, key=(lambda u: type(u).__name__))
    source = {}
    heat_utils = []
    for u in units:
        hus = u.heat_utilities
        name = u.line + ': ' + u.ID
        for i in hus: source[i] = name
        heat_utils.extend(hus)
    heat_utils.sort(key=lambda hu: hu.ID)
    
    # First table and set Type to compare with
    hu = heat_utils[0]
    Type = hu.ID
    
    # Make a list of tables, keeping all results with same Type in one table
    index = []
    table = []
    tables = []
    for hu in heat_utils:
        Type_new = hu.ID
        if Type == Type_new:
            index.append(source[hu])
            table.append((hu.duty, hu.flow, hu.cost))
        else:
            table = DataFrame(table, index=index,
                               columns=('Duty (kJ/hr)',
                                        'Flow (kmol/hr)',
                                        'Cost (USD/hr)'))
            table.columns.name = Type
            tables.append(table)
            index = [source[hu]]
            table = [(hu.duty, hu.flow, hu.cost)]
            Type = Type_new
    
    # Add the last table
    table = DataFrame(table, index=index,
                               columns=('Duty (kJ/hr)', 'Flow (kmol/hr)', 'Cost (USD/hr)'))
    table.columns.name = Type
    tables.append(table)
    return tables
    

def power_utilities_table(units):
    # power_units = units_of_measure.get('power') or 'kW'
    # cost_units = units_of_measure.get('cost') or 'USD/hr'
    # power_factor = factor('kW', power_units)
    # cost_factor =  factor('USD/hr', cost_units)
    # Sort power utilities by unit type
    units = sorted(units, key=(lambda u: type(u).__name__))
    power_utilities = [u.power_utility for u in units if u.power_utility]
    array = np.zeros([len(power_utilities), 2])
    for i, pu in enumerate(power_utilities):
        array[i] = (pu.power, pu.cost)
    return DataFrame(array, index=[f'{type(u).__name__}: {u.ID}' for u in units if u.power_utility],
                     columns=(f'Power (kW)', f'Cost (USD/hr)'))


# %% Streams

def stream_table(streams, flow='kg/min', **props) -> 'DataFrame':
    """Return a stream table as a pandas DataFrame object.

    **Parameters**:

         **streams:** array_like[Stream]
        
         **Flow:** [str] Units for flow rate

         **props:** [str] Additional stream properties and units as key-value pairs
    
    """
    # Get correct flow attributes
    flow_dim = Q_(0, flow).dimensionality
    if flow_dim == mol_flow_dim:
        flow_attr = 'mol'
    elif flow_dim == mass_flow_dim:
        flow_attr = 'mass'
    elif flow_dim == vol_flow_dim:
        flow_attr = 'vol'
    else:
        raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
    
    # Prepare rows and columns
    streams.discard(None)
    ss = [*streams]
    ss.sort(key=lambda s: s.ID)
    species = ss[0]._IDs
    n = len(ss)
    m = len(species)
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
        sources[j] = str(s.source)
        sinks[j] = str(s.sink)
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
        flow_j = getattr(s, flow_attr)
        flows[j] = net_j = sum(flow_j)
        fracs[:,j] = flow_j/net_j if net_j > 0 else 0
        i = 0
        for attr, unit in props.items():
            prop_molar_data[i, j] = getattr(s, attr)
            i += 1
    
    # Set the right units
    units = Stream.units
    flows = Q_(flows, units[flow_attr]); flows.ito(flow); flows = flows.magnitude
    i = 0
    prop_molar_keys = p*[None]
    for attr, unit in props.items():
        p = prop_molar_data[i]
        p = Q_(p, units[attr]); p.ito(unit); p = p.magnitude
        prop_molar_keys[i] = f'{attr} ({unit})'
        i += 1
    
    # Add spaces for readability
    species = list(species).copy()
    for i in range(m):
        species[i] = '- ' + species[i]
    
    # Make data frame object
    index = ('Source', 'Sink', 'Phase')  + tuple(prop_molar_keys) + (f'flow ({flow})', 'Composition:') + tuple(species)
    return DataFrame(array, columns=IDs, index=index)


