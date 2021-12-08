# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import pandas as pd
from warnings import warn
import openpyxl
import thermosteam as tmo
import biosteam as bst
from .._heat_utility import HeatUtility
import os

DataFrame = pd.DataFrame
ExcelWriter = pd.ExcelWriter

__all__ = ('stream_table', 'cost_table', 'save_system_results',
           'save_report', 'unit_result_tables', 'heat_utility_tables',
           'power_utility_table', 'tables_to_excel', 'voc_table',
           'lca_displacement_allocation_table', 
           'lca_inventory_table',
           'lca_property_allocation_factor_table',
           'lca_displacement_allocation_factor_table',
           'FOCTableBuilder')

def _stream_key(s): # pragma: no coverage
    num = s.ID[1:]
    if num.isnumeric(): return int(num)
    else: return -1

def _reformat(name):
    name = name.replace('_', ' ')
    if name.islower(): name= name.capitalize()
    return name

# %% Detailed TEA tables

class FOCTableBuilder:
    __slots__ = ('index', 'data', 'costs')
    
    def __init__(self):
        self.index = []
        self.data = []
        self.costs = []
        
    def entry(self, index, costs, notes='-'):
        self.index.append(index)
        self.data.append([notes, *costs])
        self.costs.append(costs)
    
    def table(self, names):
        data = self.data
        index = self.index.copy()
        index.append('Fixed operating cost (FOC)')
        data.append(("", *sum(self.costs)))
        return pd.DataFrame(np.array(data), 
                            index=index,
                            columns=(
                                'Notes',
                                *[i + '\n[MM$ / yr]' for i in names],
                            )
        )

# %% Multiple system tables

def voc_table(systems, product_IDs, system_names=None):
    # Not ready for users yet
    isa = isinstance
    if isa(systems, bst.System): systems = [systems]
    other_utilities_dct = {}
    other_byproducts_dct = {}
    prices = bst.stream_utility_prices
    kg_per_ton = 907.185
    def getsubdct(dct, name):
        if name in dct:
            subdct = dct[name]
        else:
            dct[name] = subdct = {}
        return subdct
    
    for sys in systems:
        electricity_cost = sys.power_utility.cost * sys.operating_hours
        if electricity_cost > 0.: 
            dct = getsubdct(other_utilities_dct, 'Electricity')
            dct[sys] = (f"{bst.PowerUtility.price} $/kWh", electricity_cost)
        else:
            dct = getsubdct(other_byproducts_dct, 'Electricity production')
            dct[sys] = (f"{bst.PowerUtility.price} $/kWh", -electricity_cost)
        inlet_flows = sys.get_inlet_utility_flows()
        for name, flow in inlet_flows.items():
            dct = getsubdct(other_utilities_dct, name)
            price = prices[name]
            dct[sys] = (price * kg_per_ton, price * flow)
        outlet_flows = sys.get_outlet_utility_flows()
        for name, flow in outlet_flows.items():
            dct = getsubdct(other_byproducts_dct, name)
            price = prices[name]
            dct[sys] = (price * kg_per_ton, price * flow)
        
    def reformat(name):
        name = name.replace('_', ' ')
        if name.islower(): name= name.capitalize()
        return name
    
    feeds = sorted({i.ID for i in sum([i.feeds for i in systems], []) if i.price})
    coproducts = sorted({i.ID for i in sum([i.products for i in systems], []) if i.price and i.ID not in product_IDs})
    system_heat_utilities = [bst.HeatUtility.sum_by_agent(sys.heat_utilities) for sys in systems]
    other_utilities = sorted(other_utilities_dct)
    other_byproducts = sorted(other_byproducts_dct)
    heating_agents = sorted(set(sum([[i.agent.ID for i in hus if i.cost and i.flow * i.duty > 0. and abs(i.flow) > 1e-6] for hus in system_heat_utilities], [])))
    cooling_agents = sorted(set(sum([[i.agent.ID for i in hus if i.cost and i.flow * i.duty < 0. and abs(i.flow) > 1e-6] for hus in system_heat_utilities], [])))
    index = {j: i for (i, j) in enumerate(feeds + heating_agents + cooling_agents + other_utilities + coproducts + other_byproducts)}
    table_index = [*[('Raw materials', reformat(i)) for i in feeds],
                   *[('Heating utilities', reformat(i)) for i in heating_agents],
                   *[('Cooling utilities', reformat(i)) for i in cooling_agents],
                   *[('Other utilities', i) for i in other_utilities],
                   *[('By-products and credits', reformat(i)) for i in coproducts],
                   *[('By-products and credits', i) for i in other_byproducts]]
    table_index.append(('Variable operating cost', ''))
    N_cols = len(systems) + 1
    N_rows = len(table_index)
    data = np.zeros([N_rows, N_cols], dtype=object)
    N_coproducts = len(coproducts) + len(other_byproducts)
    for col, sys in enumerate(systems):
        for stream in sys.feeds + sys.products:
            if stream.ID in product_IDs: continue
            cost = sys.get_market_value(stream)
            if cost:
                ind = index[stream.ID]
                data[ind, 0] = stream.price * kg_per_ton # USD / ton
                data[ind, col + 1] = cost / 1e6 # million USD / yr
        for hu in system_heat_utilities[col]:
            try:
                ind = index[hu.agent.ID]
            except:
                continue
            cost = sys.operating_hours * hu.cost
            price = ""
            if hu.agent.heat_transfer_price:
                price += f"{hu.agent.heat_transfer_price} USD/kJ"
            if hu.agent.regeneration_price:
                if price: price += ", "
                price += f"{hu.agent.regeneration_price} USD/kmol"
            data[ind, 0] = price 
            data[ind, col + 1] = cost / 1e6 # million USD / yr
        for sysdct in (other_byproducts_dct, other_utilities_dct):
            for i, dct in sysdct.items():
                if sys not in dct: continue
                price, cost = dct[sys]
                ind = index[i]
                data[ind, 0] = price
                data[ind, col + 1] = cost / 1e6 # million USD / yr
    N_consumed = N_rows - N_coproducts - 1
    data[-1, 1:] = data[:N_consumed, 1:].sum(axis=0) - data[N_consumed:, 1:].sum(axis=0)
    if system_names is None:
        system_names = [i.ID for i in systems]
    columns = [i + " [MM$/yr]" for i in system_names]
    return pd.DataFrame(data, 
                        index=pd.MultiIndex.from_tuples(table_index),
                        columns=('Price [$/ton]', *columns))

def lca_inventory_table(systems, key, items=(), system_names=None):
    items = frozenset(items)
    isa = isinstance
    if isa(systems, bst.System): systems = [systems]
    PowerUtility = bst.PowerUtility
    other_utilities = []
    other_byproducts = []
    process_items = []
    other_values = {}
    def set_value(name, sys, value):
        if name in other_values:
            dct = other_values[name]
        else:
            other_values[name] = dct = {}
        dct[sys] = value
        
    for sys in systems:
        electricity_consumption = sys.power_utility.rate * sys.operating_hours
        if electricity_consumption > 0.: 
            if 'Electricity [kWhr/yr]' not in other_utilities: other_utilities.append('Electricity [kWhr/yr]')
            cf = PowerUtility.get_CF(key, production=False)
            if cf:
                set_value('Electricity [kWhr/yr]', sys, electricity_consumption)
        else:
            if 'Electricity [kWhr/yr]' not in other_byproducts: other_byproducts.append('Electricity [kWhr/yr]')
            cf = PowerUtility.get_CF(key, consumption=False)
            if cf:
                set_value('Electricity [kWhr/yr]', sys, -electricity_consumption)
        try: process_impact_items = sys.process_impact_items[key]
        except: continue
        for item in process_impact_items:
            if item.name not in process_items: process_items.append(item.name)
            value = item.inventory()
            if item.basis != 'kg':
                value = f"{value} [{item.basis}/yr]"
            set_value(item.name, sys, value)
        
    def reformat(name):
        name = name.replace('_', ' ')
        if name.islower(): name= name.capitalize()
        return name
    
    feeds = sorted({i.ID for i in sum([i.feeds for i in systems], []) if key in i.characterization_factors})
    coproducts = sorted({i.ID for i in sum([i.products for i in systems], []) if key in i.characterization_factors or i in items})
    system_heat_utilities = [bst.HeatUtility.sum_by_agent(sys.heat_utilities) for sys in systems]
    input_heating_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty > 0. and i.flow > 1e-6] for hus in system_heat_utilities], [])))
    input_cooling_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty < 0. and i.flow > 1e-6] for hus in system_heat_utilities], [])))
    output_heating_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty > 0. and i.flow < -1e-6] for hus in system_heat_utilities], [])))
    output_cooling_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty < 0. and i.flow < -1e-6] for hus in system_heat_utilities], [])))
    index = {j: i for (i, j) in enumerate(feeds + input_heating_agents + input_cooling_agents + other_utilities + coproducts + other_byproducts + output_heating_agents + output_cooling_agents + process_items)}
    table_index = [*[('Inputs', reformat(i)) for i in feeds + input_heating_agents + input_cooling_agents + other_utilities],
                   *[('Outputs', reformat(i)) for i in coproducts + other_byproducts + output_heating_agents + output_cooling_agents],
                   *[(i, '') for i in process_items]]
    N_cols = len(systems)
    N_rows = len(table_index)
    data = np.zeros([N_rows, N_cols], dtype=object)
    for col, sys in enumerate(systems):
        for stream in sys.feeds + sys.products:
            try: ind = index[stream.ID]
            except: continue
            data[ind, col] = sys.get_mass_flow(stream)
        for hu in system_heat_utilities[col]:
            try: ind = index[hu.agent.ID]
            except: continue
            flow, units = hu.get_inventory()
            if flow:
                flow = sys.operating_hours * flow
                data[ind, col] = f"{flow} [{units}]"
        for i, subdct in other_values.items():
            if sys not in subdct: continue
            value = subdct[sys]
            ind = index[i]
            data[ind, col] = value
    if system_names is None and len(systems) == 1:
        columns = ["Inventory [kg/yr]"]
    else:
        if system_names is None:
            system_names = [i.ID for i in systems]
        sys_units = " inventory [kg/yr]"
        columns = [i + sys_units for i in system_names]
    return pd.DataFrame(data, 
                        index=pd.MultiIndex.from_tuples(table_index),
                        columns=columns)

def lca_displacement_allocation_table(systems, key, items, 
                                      item_name=None, system_names=None):
    # Not ready for users yet
    isa = isinstance
    if isa(systems, bst.System): systems = [systems]
    PowerUtility = bst.PowerUtility
    other_utilities = []
    other_byproducts = []
    process_inventory = []
    other_values = {}
    impact_units = str(bst.settings.get_impact_indicator_units(key))
    def set_value(name, sys, CF, value):
        if name in other_values:
            dct = other_values[name]
        else:
            other_values[name] = dct = {}
        dct[sys] = (CF, value)
        
    for sys in systems:
        electricity_consumption = sys.power_utility.rate * sys.operating_hours
        if electricity_consumption > 0.: 
            if 'Electricity' not in other_utilities: other_utilities.append('Electricity')
            cf = PowerUtility.get_CF(key, production=False)
            set_value('Electricity', sys, f"{cf} {impact_units}/kWhr", electricity_consumption * cf)
        else:
            if 'Electricity' not in other_byproducts: other_byproducts.append('Electricity')
            cf = PowerUtility.get_CF(key, consumption=False)
            set_value('Electricity', sys, f"{cf} {impact_units}/kWhr", - electricity_consumption * cf)
        try: process_impact_items = sys.process_impact_items[key]
        except: continue
        for item in process_impact_items:
            if item.name not in process_inventory: process_inventory.append(item.name)
            value = item.impact()
            CF = item.CF
            basis = item.basis
            if basis != 'kg':
                CF = f"{CF} [{impact_units}/{basis}"
                value = f"{value} [{basis}/yr]"
            set_value(item.name, sys, CF, value)
    
    if item_name is None: item_name = items[0].ID.replace('_', ' ')
    feeds = sorted({i.ID for i in sum([i.feeds for i in systems], []) if key in i.characterization_factors})
    coproducts = sorted({i.ID for i in sum([i.products for i in systems], []) if key in i.characterization_factors})
    system_heat_utilities = [bst.HeatUtility.sum_by_agent(sys.heat_utilities) for sys in systems]
    input_heating_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty > 0. and i.flow > 1e-6] for hus in system_heat_utilities], [])))
    input_cooling_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty < 0. and i.flow > 1e-6] for hus in system_heat_utilities], [])))
    output_heating_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty > 0. and i.flow < -1e-6] for hus in system_heat_utilities], [])))
    output_cooling_agents = sorted(set(sum([[i.agent.ID for i in hus if (key, i.agent.ID) in i.characterization_factors and i.flow * i.duty < 0. and i.flow < -1e-6] for hus in system_heat_utilities], [])))
    inputs = [*feeds, *input_heating_agents, *input_cooling_agents, *other_utilities]
    outputs = [*coproducts, *other_byproducts, *output_heating_agents, *output_cooling_agents]
    keys = [*inputs, 'Total inputs', *outputs, 'Total outputs displaced', *process_inventory, 'Total']
    index = {j: i for (i, j) in enumerate(keys)}
    table_index = [*[('Inputs', _reformat(i)) for i in feeds + input_heating_agents + input_cooling_agents + other_utilities],
                   ('Total inputs', ''),
                   *[('Outputs displaced', _reformat(i)) for i in coproducts + other_byproducts + output_heating_agents + output_cooling_agents],
                   ('Total outputs displaced', ''),
                   *[("Process impacts", i) for i in process_inventory],
                   ("Total", '')]
    N_cols = len(systems) + 1
    N_rows = len(table_index)
    data = np.zeros([N_rows, N_cols], dtype=object)
    for col, sys in enumerate(systems):
        item_flow = sys.get_mass_flow(items[col])
        for stream in sys.feeds + sys.products:
            try: ind = index[stream.ID]
            except: continue
            data[ind, 0] = stream.characterization_factors[key]
            impact = sys.get_material_impact(stream, key)
            data[ind, col + 1] = impact / item_flow
        for hu in system_heat_utilities[col]:
            try: ind = index[hu.agent.ID]
            except: continue
            impact = sys.operating_hours * hu.get_impact(key)
            cf, basis = hu.characterization_factors[hu.agent.ID, key]
            data[ind, 0] = f"{cf} [{impact_units}{basis}]"
            data[ind, col + 1] = impact / item_flow
        for i, subdct in other_values.items():
            if sys not in subdct: continue
            cf, impact = subdct[sys]
            ind = index[i]
            data[ind, 0] = cf
            data[ind, col + 1] = impact / item_flow
    total_inputs_index = index['Total inputs']
    data[total_inputs_index, 1:] = total_inputs = data[:total_inputs_index, 1:].sum(axis=0)
    data[total_inputs_index, 0] = ''
    total_outputs_index = index['Total outputs displaced']
    data[total_outputs_index, 1:] = total_outputs = data[total_inputs_index + 1:total_outputs_index, 1:].sum(axis=0)
    data[total_outputs_index, 0] = ''
    data[-1, 1:] = total_inputs - total_outputs + data[total_outputs_index + 1:-1, 1:].sum(axis=0)
    data[-1, 0] = ''
    if system_names is None and len(systems) == 1:
        columns = [f"{key} [{impact_units}/kg*{item_name}]"]
    else:
        if system_names is None:
            system_names = [i.ID for i in systems]
        sys_units = f" {key} [{impact_units}/kg*{item_name}]"
        columns = [i + sys_units for i in system_names]
    return pd.DataFrame(data, 
                        index=pd.MultiIndex.from_tuples(table_index),
                        columns=(f'Characterization factor [{impact_units}/kg]', *columns))

def lca_property_allocation_factor_table(
        systems, property, units, system_names=None
    ):
    system_allocation_factors = [i.get_property_allocation_factors(property, units) for i in systems]
    table_index = sorted(set(sum([tuple(i) for i in system_allocation_factors], ())))
    index = {j: i for i, j in enumerate(table_index)}
    N_cols = len(systems)
    N_rows = len(table_index)
    data = np.zeros([N_rows, N_cols], dtype=object)
    for col, factors in enumerate(system_allocation_factors):
        for key, value in factors.items():
            data[index[key], col] = value
    if system_names is None and len(systems) == 1:
        columns = [f"{_reformat(property)} allocation factors"]
    else:
        if system_names is None:
            system_names = [i.ID for i in systems]
        sys_units = f" {property.replace('_', ' ')} allocation factors"
        columns = [i + sys_units for i in system_names]
    return pd.DataFrame(data, 
                        index=table_index,
                        columns=columns)

def lca_displacement_allocation_factor_table(
        systems, items, key, system_names=None
    ):
    system_allocation_factors = [i.get_displacement_allocation_factors(j, key) for i, j in zip(systems, items)]
    table_index = sorted(set(sum([tuple(i) for i in system_allocation_factors], ())))
    index = {j: i for i, j in enumerate(table_index)}
    N_cols = len(systems)
    N_rows = len(table_index)
    data = np.zeros([N_rows, N_cols], dtype=object)
    for col, factors in enumerate(system_allocation_factors):
        for key, value in factors.items():
            data[index[key], col] = value
    if system_names is None and len(systems) == 1:
        columns = ["Displacement allocation factors"]
    else:
        if system_names is None:
            system_names = [i.ID for i in systems]
        sys_units = " displacement allocation factors"
        columns = [i + sys_units for i in system_names]
    return pd.DataFrame(data, 
                        index=table_index,
                        columns=columns)

def lca_allocation_method_results(
        systems, key, property_args, items, item_name 
    ):
    pass

# %% Helpful functions

def tables_to_excel(tables, writer, sheet='Sheet1', n_row=1, row_spacing=2): 
    """
    Save a list of tables as an excel file and return the row number at which
    another consecutive table would start.
    
    Parameters
    ----------
    tables : iterable[pandas.DataFrame]
        Tables to save to excel.
    writer : pandas.ExcelWritter
        Writer to manage data stream to excel.
    sheet : str
        Name of sheet to save data.
    n_row : int
        Row number to begin saving data.
    row_spacing : int
        Number of rows between tables.
    
    Returns
    -------
    n_row : int
        Row number for next table.
        
    """
    row_spacing += 1 # Account for Python index offset
    for t in tables:
        label = t.columns.name
        t.to_excel(writer, sheet, 
                   startrow=n_row, index_label=label)
        n_row += len(t.index) + row_spacing
    return n_row

# %% Units

def save_report(system, file='report.xlsx', dpi='300', tea=None, **stream_properties): 
    """
    Save a system report as an xlsx file.
    
    Parameters
    ----------
    file : str
        File name to save report
    dpi : str, optional
        Resolution of the flowsheet. Defaults to '300'
    tea : TEA, optional
        Object for techno-economic analysis and cashflows. Defaults to the
        TEA object linked to the system.
    **stream_properties : str
        Additional stream properties and units as key-value pairs (e.g. T='degC', flow='gpm', H='kW', etc..)
        
    """
    writer = ExcelWriter(file)
    units = sorted(system.units, key=lambda x: x.line)
    cost_units = [i for i in units if i._design or i._cost]
    try:
        system.diagram('thorough', file='flowsheet', dpi=str(dpi), format='png')
    except:
        diagram_completed = False
        warn(RuntimeWarning('failed to generate diagram through graphviz'), stacklevel=2)
    else:
        import PIL.Image
        try:
            # Assume openpyxl is used
            worksheet = writer.book.create_sheet('Flowsheet')
            flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
            worksheet.add_image(flowsheet, anchor='A1')
        except PIL.Image.DecompressionBombError:
            PIL.Image.MAX_IMAGE_PIXELS = int(1e9)
            flowsheet = openpyxl.drawing.image.Image('flowsheet.png')
            worksheet.add_image(flowsheet, anchor='A1')
        except:
            # Assume xlsx writer is used
            try:
                worksheet = writer.book.add_worksheet('Flowsheet')
            except:
                breakpoint()
            worksheet.insert_image('A1', 'flowsheet.png')
        diagram_completed = True
    
    if tea is None: tea = system.TEA
    if tea:
        tea = system.TEA
        cost = cost_table(tea)
        cost.to_excel(writer, 'Itemized costs')
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
    for chemicals, streams in streams_by_chemicals.items():
        stream_tables.append(stream_table(streams, chemicals=chemicals, T='K', **stream_properties))
    tables_to_excel(stream_tables, writer, 'Stream table')
    
    # Heat utility tables
    heat_utilities = heat_utility_tables(cost_units)
    n_row = tables_to_excel(heat_utilities, writer, 'Utilities')
    
    # Power utility table
    power_utility = power_utility_table(cost_units)
    power_utility.to_excel(writer, 'Utilities', 
                           index_label='Electricity',
                           startrow=n_row)
    
    # General desing requirements
    results = unit_result_tables(cost_units)
    tables_to_excel(results, writer, 'Design requirements')
    
    # Reaction tables
    reactions = unit_reaction_tables(units)
    tables_to_excel(reactions, writer, 'Reactions')
    
    writer.save()
    if diagram_completed: os.remove("flowsheet.png")

save_system_results = save_report

def unit_reaction_tables(units):
    isa = isinstance
    tables = []
    rxntypes = (tmo.Reaction, tmo.ReactionSet)
    for u in units:
        all_reactions = {rxn for rxn in u.__dict__.values() if isa(rxn, rxntypes)}
        for rxn in tuple(all_reactions):
            if hasattr(rxn, '_parent'):
                if rxn._parent in all_reactions: all_reactions.discard(rxn)
            elif hasattr(rxn, '_parent_index'):
                parent, index = rxn._parent_index
                if parent in all_reactions: all_reactions.discard(rxn)
        for rxn in all_reactions:
            df = rxn.to_df()
            df.columns.name = '-'.join([u.ID, u.line])
            tables.append(df)
    return tables

def unit_result_tables(units,
                       include_utilities=False, 
                       include_total_cost=False,
                       include_installed_cost=False): 
    """
    Return a list of results tables for each unit type.

    Parameters
    ----------
    units : iterable[Unit]
        
    Returns
    -------
    tables : list[DataFrame]
    
    """
    units = sorted(units, key=(lambda u: u.line))
    
    # Organize units by units of measure:
    organized = {}
    for u in units:
        uom = (*u._units.keys(), u.line)
        if uom in organized: organized[uom].append(u)
        else: organized[uom] = [u]
    
    # Make a list of tables, keeping all results with same keys in one table
    tables = []
    key = lambda u: u.ID
    for all_units in organized.values():
        # First table with units of measures
        all_units = sorted(all_units, key=key)
        u, *units = all_units
        empty_heat_utilities = []
        key_hook = None
        if all([i.line == 'Heat Exchanger' for i in all_units]):
            key_hook = lambda key: ('Purchase cost', 'Heat Exchanger') if key[0]=='Purchase cost' else key
            all_heat_utilities = sum([i.heat_utilities for i in all_units], ())
            for i in set([i.agent for i in all_heat_utilities]):
                if not i: continue
                heat_utility = HeatUtility()
                heat_utility.load_agent(i)
                empty_heat_utilities.append(heat_utility)
        empty_heat_utilities = tuple(empty_heat_utilities)
        table = u.results(include_utilities=include_utilities,
                          include_total_cost=include_total_cost,
                          include_installed_cost=include_installed_cost,
                          include_zeros=False,
                          external_utilities=empty_heat_utilities,
                          key_hook=key_hook)
        if table is None: continue
        for u in units:
            table[u.ID] = u.results(with_units=False, 
                                    include_utilities=include_utilities,
                                    include_total_cost=include_total_cost,
                                    include_installed_cost=include_installed_cost,
                                    include_zeros=False,
                                    external_utilities=empty_heat_utilities,
                                    key_hook=key_hook)
        table.columns.name = (u.line, '')
        tables.append(table)
    return tables
    
def cost_table(tea): 
    """
    Return a cost table as a pandas DataFrame object.

    Parameters
    ----------
    units : iterable[Unit]
        
    Returns
    -------
    table : DataFrame

    """
    columns = ('Unit operation',
               'Purchase cost (10^6 USD)',
               'Utility cost (10^6 USD/yr)')
    units = sorted([i for i in tea.system.units if i._design or i._cost], key=lambda x: x.line)
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
        df['Installed cost (10^6 USD)'] = [u.installed_cost / 1e6 for u in units]
    
    return df

def heat_utility_tables(units): 
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
    

def power_utility_table(units): 
    """
    Return a pandas DataFrame object of power utilities.
    
    Parameters
    ----------
    units : iterable[Unit]
        
    """
    # Sort power utilities by unit type
    units = sorted(units, key=(lambda u: type(u).__name__))
    units = [u for u in units if u.power_utility]
    power_utilities = [u.power_utility for u in units]
    length = len(power_utilities)
    data = []
    for i, u, pu in zip(range(length), units, power_utilities):
        data.append((u.line, pu.rate, pu.cost))
    return DataFrame(data, index=[u.ID for u in units if u.power_utility],
                     columns=('Unit Operation', 'Rate (kW)', 'Cost (USD/hr)'))

# def lca_tables(sys):
#     all_feeds = [i for i in sys.feeds if i.characterization_factors]
#     all_products = [i for i in sys.products if i.characterization_factors]
#     net_electricity = bst.PowerUtility.sum([i.power_utility for i in sys.cost_units])
#     keys = list(net_electricity.characterization_factors)
#     keys = set(sum([list(i.characterization_factors) for i in (all_feeds + all_products)], keys))
#     tables = []
#     for key in keys:
#         feeds = [i for i in all_feeds if key in i.characterization_factors]
#         mass_flows = [i.F_mass for i in feeds]
        
#         df = DataFrame(data, index=[])
#         tables.append(df)

# %% Streams

def stream_table(streams, flow='kg/hr', percent=True, chemicals=None, **props):
    """
    Return a stream table as a pandas DataFrame object.

    Parameters
    ----------
    streams : array_like[Stream]
    flow : str
        Units for flow rate.
    props : str
        Additional stream properties and units as key-value pairs
    
    """
    
    # Prepare rows and columns
    ss = sorted(sorted([i for i in streams if i.ID], key=lambda i: i.ID), key=_stream_key)
    if not chemicals: 
        all_chemicals = tuple(set([i.chemicals for i in ss]))
        sizes = [(i, chemical.size) for i, chemical in enumerate(all_chemicals)]
        index, size = max(sizes, key=lambda x: x[1])
        chemicals = all_chemicals[index]
    n = len(ss)
    m = chemicals.size
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
        sources[j] = s.source.ID if s.source else '-'
        sinks[j] = s.sink.ID if s.sink else '-'
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
        flow_j = s.get_flow(units=flow)
        flows[j] = net_j = flow_j.sum()
        if percent: net_j /= 100.
        fracs_j = flow_j/net_j if net_j > 1e-24 else 0
        if s.chemicals is chemicals:
            fracs[:, j] = fracs_j
        else:
            fracs[:, j] = 0.
            fracs[chemicals.get_index(s.chemicals.IDs), j] = fracs_j
        i = 0
        for attr, units in props.items():
            prop_molar_data[i, j] = s.get_property(attr, units)
            i += 1
    index = (
        'Source', 
        'Sink',
        'Phase', 
        *[f'{attr} ({units})' for attr, units in props.items()], 
        f'flow ({flow})',
        ('Composition [%]:' if percent else 'Composition:'),
        *chemicals.IDs
    )
    return DataFrame(array, columns=IDs, index=index)


