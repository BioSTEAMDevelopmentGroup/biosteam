# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 09:48:34 2018

@author: yoelr
"""
from biosteam.stream import Stream, mol_flow_dim, mass_flow_dim, vol_flow_dim
from biosteam import Q_, pd, np
from biosteam.exceptions import DimensionError

DataFrame = pd.DataFrame
ExcelWriter = pd.ExcelWriter

# %% Helpful functions

def units_wt_cost(units):
    """Remove units that do not have capital or operating cost."""
    r = 0
    for i in range(len(units)):
        i -= r
        unit = units[i]
        ci = unit.capital_cost
        co = unit.operating_cost
        if not (ci or co):
            del units[i]
            r += 1

def unit_sort(units):
    """Sort by first grouping units by type and then order groups by maximum capital."""
    unit_types = list({type(u).__name__ for u in units})
    type_dict = {ut:[] for ut in unit_types}
    for u in units:
        ut = type(u).__name__
        type_dict[ut].append(u)
    units = []
    unit_types.sort(key=lambda ut: max(u.capital_cost for u in type_dict[ut]), reverse=True)
    for key in unit_types:
        ulist = type_dict[key]
        ulist.sort(key=lambda x: x.capital_cost, reverse=True)
        units += ulist
    return units


# %% Units

def unit_report(units, name='Unit Report.xlsx'):
    """Return a unit report as a pandas DataFrame object.

    **Parameters**:

        **units:** array_like[Unit]
         
    **Returns**
    
        **unit_report:** [DataFrame] Units as columns with all results as rows.

    """
    units_wt_cost(units)
    units = unit_sort(units)
    ur = units[0].results
    keys = tuple(ur.nested_keys())
    report = ur.table()
    writer = ExcelWriter(name)
    del units[0]
    n_row = 1
    for u in units:
        ur = u.results
        new_keys = tuple(ur.nested_keys())
        if new_keys == keys:
            report = pd.concat( (report, ur.table(with_units=False)), axis=1)
        else:
            label = report.columns.name
            report.to_excel(writer, startrow=n_row, index_label=label)
            n_row += len(report.index) + 2
            report = ur.table()
            keys = new_keys
    
    
def unit_table(units) -> 'DataFrame':
    """Return a unit table as a pandas DataFrame object.

    **Parameters**:

        **units:** array_like[Unit]
         
    **Returns**
    
        **unit_table:** [DataFrame] Units as indexes with the following columns
            * 'Unit Type': Type of unit
            * 'CI (10^6 USD)': Capital investment
            * 'OC (10^6 USD/yr)': Annual operating cost
            * '% TCI': Percent of total capital investment
            * '% TOC': Percent of operating cost

    """
    columns = ('Type',
              f'CI (10^6 USD)',
              f'OC (10^6 USD/yr)',
              f'% TCI',
              f'% TOC')
    units_wt_cost(units)
    units = unit_sort(units)
    
    # Initialize data
    N_units = len(units)
    array = np.empty((N_units, 5), dtype=object)
    IDs = []
    types = array[0:, 0]
    C_cap = array[0:, 1]
    C_op = array[0:, 2]
    
    # Get data
    for i in range(N_units):
        unit = units[i]
        types[i] = type(unit).__name__
        C_cap[i] = unit.capital_cost * unit.lang_factor / 1e6
        C_op[i] = unit.operating_cost * unit.operating_days * 24  / 1e6
        IDs.append(unit.ID)
    
    # Make data frame
    array[0:, 3] = C_cap/sum(C_cap)*100
    array[0:, 4] = C_op/sum(C_op)*100
    return DataFrame(array, columns=columns, index=IDs)
    
# %% Streams

def stream_table(streams, Flow='kmol/hr', **props) -> 'DataFrame':
    """Return a stream table as a pandas DataFrame object.

    **Parameters**:

         **streams:** array_like[Stream]
        
         **Flow:** [str] Units for flow rate

         **props:** [str] Additional stream properties and units as key-value pairs
    
    """
    # Get correct flow attributes
    flow_dim = Q_(0, Flow).dimensionality
    if flow_dim == mol_flow_dim:
        flow_attr = 'mol'
    elif flow_dim == mass_flow_dim:
        flow_attr = 'mass'
    elif flow_dim == vol_flow_dim:
        flow_attr = 'vol'
    else:
        raise DimensionError(f"Dimensions for flow units must be in molar, mass or volumetric flow rates, not '{flow_dim}'.")
    
    # Prepare rows and columns
    ss = streams
    species = type(ss[0])._specie_IDs
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
        sources[j] = s.source[0]
        sinks[j] = s.sink[0]
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
    flows = Q_(flows, units[flow_attr]); flows.ito(Flow); flows = flows.magnitude
    i = 0
    prop_molar_keys = p*[None]
    for attr, unit in props.items():
        p = prop_molar_data[i]
        p = Q_(p, units[attr]); p.ito(unit); p = p.magnitude
        prop_molar_keys[i] = f'{attr} ({unit})'
        i += 1
    
    # Add spaces for readability
    species = species.copy()
    for i in range(m):
        species[i] = '- ' + species[i]
    
    # Make data frame object
    index = ('Source', 'Sink', 'Phase')  + tuple(prop_molar_keys) + (f'Flow ({Flow})', 'Composition:') + tuple(species)
    return DataFrame(array, columns=IDs, index=index)