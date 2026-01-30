# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
import numpy as np
import pandas as pd
from .._heat_utility import HeatUtility

__all__ = (
    'group_units_by_available_chemicals',
    'is_storage_unit',
    'rename_unit',
    'rename_units',
    'group_by_area',
    'heat_exchanger_utilities_from_units',
    'units_with_costs',
    'units_with_heat_exchangers',
    'filter_out_heat_utility_savings',
    'get_utility_flow',
    'get_utility_duty',
    'get_power_utilities',
    'get_heat_utilities',
    'get_purchase_cost',
    'get_installed_cost',
    'get_cooling_duty',
    'get_heating_duty',
    'get_electricity_consumption',
    'get_electricity_production',
    'get_tea_results',
    'group_by_lines',
    'group_by_types',
    'filter_by_types',
    'filter_by_lines',
    'volume_of_chemical_in_units',
    'set_construction_material',
    'set_construction_material_to_stainless_steel',
    'set_construction_material_to_carbon_steel',
    'get_OSBL',
    'heat_exchanger_operation',
    'connect_by_ID',
    'default_utilities',
    'default_CEPCI',
    'default',
)

def connect_by_ID(units):
    """Connect unit inlets and outlets where streams have the same ID."""
    connections = []
    inlets = {}
    outlets = {}
    units = tuple(units)
    for u in units:
        for i, s in enumerate(u._ins):
            ID = s.ID
            if ID in outlets: connections.append(ID)
            inlets[ID] = (u, i)
        for i, s in enumerate(u._outs):
            ID = s.ID
            if ID in inlets: connections.append(ID)
            outlets[ID] = (u, i)
    for ID in connections:
        upstream, index = outlets[ID]
        stream = upstream.outs[index]
        downstream, index = inlets[ID]
        downstream.ins[index] = stream

def group_units_by_available_chemicals(units):
    groups = {}
    for u in units:
        available_chemicals = tuple(u.get_available_chemicals())
        if available_chemicals in groups:
            groups[available_chemicals].append(u)
        else:
            groups[available_chemicals] = [u]
    return groups
        
def is_storage_unit(unit):
    return (
        ('storage' in unit.line.lower() 
         or isinstance(unit, bst.StorageTank)
         or 'storage' in unit.__class__.__name__.lower()) 
        and (unit.ins[0].isfeed() or unit.outs[0].isproduct())
    )

def area_convention_number(unit_registry, letter, number):
    ID = letter + str(number)
    if ID in unit_registry:
        number += 1
        return area_convention_number(unit_registry, letter, number)
    else:
        return number

def register_by_area_convention(unit, area, unit_registry):
    letter = unit.ticket_name
    number = area_convention_number(unit_registry, letter, area[letter])
    area[letter] = number + 1 
    ID = letter + str(number)
    unit_registry.register(ID, unit)

def rename_unit(unit, area):
    """
    Rename unit according to area convention.
    
    Parameters
    ----------
    unit : :class:`biosteam.Unit`
        Unit to rename.
    area : int
        Ending number.
        
    Notes
    -----
    The area convention follows "{letter}{area + number}" where the letter depends on
    the unit operation as follows:
    
    * C = Centrifuge
    * D = Distillation column
    * E = Evaporator
    * F = Flash tank
    * H = Heat exchange
    * K = Compressor
    * Ʞ = Turbine
    * M = Mixer
    * P = Pump (including conveying belt)
    * R = Reactor
    * S = Splitter (including solid/liquid separator)
    * T = Tank or bin for storage
    * U = Other units
    * V = Valve
    * J = Junction, not a physical unit (serves to adjust streams)
    * PS = Process specification, not a physical unit (serves to adjust streams)
    
    Examples
    --------
    >>> from biosteam import *
    >>> main_flowsheet.clear() # Remove any previous data
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> M1 = Mixer('M1')
    >>> rename_unit(M1, 200)
    >>> print(M1)
    M201
    
    """
    unit_registry = bst.main_flowsheet.unit
    unit_registry.discard(unit)
    letter = unit.ticket_name
    number = area_convention_number(unit_registry, letter, int(area) + 1) 
    ID = letter + str(number)
    unit_registry.discard(unit)
    unit_registry.register(ID, unit)

def rename_units(units, area):
    """
    Rename units according to area convention.
    
    Parameters
    ----------
    units : Set[:class:`biosteam.Unit`]
        Units to rename.
    area : int
        Ending number.
        
    Notes
    -----
    The area convention follows "{letter}{area + number}" where the letter depends on
    the unit operation as follows:
    
    * C = Centrifuge
    * D = Distillation column
    * E = Evaporator
    * F = Flash tank
    * H = Heat exchange
    * K = Compressor
    * M = Mixer
    * P = Pump (including conveying belt)
    * R = Reactor
    * S = Splitter (including solid/liquid separator)
    * T = Tank or bin for storage
    * U = Other units
    * J = Junction, not a physical unit (serves to adjust streams)
    * PS = Process specification, not a physical unit (serves to adjust streams)
    
    Examples
    --------
    >>> from biosteam import *
    >>> main_flowsheet.clear() # Remove any previous data
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> units = [Mixer(),
    ...          Mixer(),
    ...          ShortcutColumn(LHK=['Ethanol', 'Water'], k=1.2, Lr=0.8, Hr=0.9),
    ...          Flash(P=101325, V=0.5),
    ...          Pump(),
    ...          Splitter(split=0.5),
    ...          MixTank(),
    ...          StorageTank(),
    ...          MultiEffectEvaporator(P=[101325, 9e4], V=0.5)]
    >>> rename_units(units, 200)
    >>> units
    [<Mixer: M201>, <Mixer: M202>, <ShortcutColumn: D201>, <Flash: F201>, <Pump: P201>, <Splitter: S201>, <MixTank: T201>, <StorageTank: T202>, <MultiEffectEvaporator: E201>]
    
    >>> # ID conflicts are taken care of internally
    >>> mixer, *other_units = units
    >>> rename_units(other_units, 200)
    >>> units
    [<Mixer: M201>, <Mixer: M202>, <ShortcutColumn: D201>, <Flash: F201>, <Pump: P201>, <Splitter: S201>, <MixTank: T201>, <StorageTank: T202>, <MultiEffectEvaporator: E201>]
    
    """
    area = int(area)
    area += 1
    units = tuple(units)
    area_dct = {i.ticket_name: area for i in units}
    unit_registry = bst.main_flowsheet.unit
    for i in units: unit_registry.discard(i)
    for i in units: register_by_area_convention(i, area_dct, unit_registry)

def units_with_costs(units):
    """Filter units that have a cost or design method."""
    return [i for i in units if i._cost or i._design]

def units_with_heat_exchangers(units):
    """Return a list of units with heat exchangers."""
    return [u for u in units if any([i.heat_exchanger for i in u.heat_utilities])]

def heat_exchanger_utilities_from_units(units):
    """Return a list of heat utilities from all heat exchangers,
    including the condensers and boilers of distillation columns and
    flash vessel heat exchangers."""
    heat_utilities = sum([i.heat_utilities for i in units], [])
    return [i for i in heat_utilities if i.hxn_ok and i.flow > 0.]

def ID_number(ID):
    """
    Return number if ID follows naming convention "{letter}{number}" where 
    the area is an integer. Returns None if no area is found
    
    """
    for i, letter in enumerate(ID):
        if letter.isdigit(): break
    else:
        return ''
    for j, letter in enumerate(ID[i:], start=i+1):
        if not letter.isdigit(): 
            j -= 1
            break
    return ID[i:j]
    
def ID_area(ID):
    """
    Return area if ID follows naming convention "{letter}{area}{digit}{digit}" 
    where the area is an integer. Returns None if no area is found
    
    """
    number = ID_number(ID)
    N_digits = len(number)
    if N_digits < 3: return 0
    return int(number[:-2] + '00')

def group_by_area(units):
    """Create a dictionary containing lists of UnitGroup objects by area."""
    areas = [ID_area(i.ID) for i in units]
    groups = {i: [] for i in sorted(set(areas))}
    for a, u in zip(areas, units): groups[a].append(u)
    return groups

def group_by_lines(units): 
    """Return a dictionary of lists of units grouped by line."""
    groups = {i.line: [] for i in units}
    for i in units: groups[i.line].append(i)
    return groups

def group_by_types(units): 
    """Return a dictionary of lists of units grouped by type."""
    groups = {i.__class__: [] for i in units}
    for i in units: groups[i.__class__].append(i)
    return groups

def filter_by_types(units, types):
    """Filter units by type(s)."""
    isa = isinstance
    return [i for i in units if isa(i, types)]

def filter_by_lines(units, lines):
    """Filter units by line(s)."""
    return [i for i in units if i.line in lines]

def filter_out_heat_utility_savings(heat_utilities):
    return [i for i in heat_utilities if i.flow > 0.]

def get_utility_flow(heat_utilities, agent):
    """Return the total utility duty of heat utilities for given agent in GJ/hr"""
    if isinstance(agent, str): agent = HeatUtility.get_agent(agent)
    return sum([i.flow * i.agent.MW for i in heat_utilities if i.agent is agent]) / 1e3

def get_utility_duty(heat_utilities, agent):
    """Return the total utility duty of heat utilities for given agent in GJ/hr"""
    if isinstance(agent, str): agent = HeatUtility.get_agent(agent)
    return sum([i.duty for i in heat_utilities if i.agent is agent]) / 1e6 

def get_power_utilities(units):
    """Return a list of all PowerUtility objects."""
    return [i.power_utility for i in units if i.power_utility]

def get_heat_utilities(units):
    """Return a list of all HeatUtility objects."""
    return sum([i.heat_utilities for i in units], [])

def get_purchase_cost(units):
    """Return the total equipment purchase cost of all units in million USD."""
    return sum([i.purchase_cost for i in units]) / 1e6 # millions USD

def get_installed_cost(units):
    """Return the total installed equipment cost of all units in million USD."""
    return sum([i.installed_cost for i in units]) / 1e6 # millions USD

def get_cooling_duty(heat_utilities, filter_savings=True):
    """Return the total cooling duty of all heat utilities in GJ/hr."""
    return - sum([i.duty for i in heat_utilities if i.flow * i.duty < 0]) / 1e6 # GJ/hr
               
def get_heating_duty(heat_utilities, filter_savings=True):
    """Return the total heating duty of all heat utilities in GJ/hr."""
    return sum([i.duty for i in heat_utilities if i.flow * i.duty > 0]) / 1e6 # GJ/hr
          
def get_electricity_consumption(power_utilities):
    """Return the total electricity consumption of all PowerUtility objects in MW."""
    return sum([i.consumption for i in power_utilities]) / 1000 # MW

def get_electricity_production(power_utilities):
    """Return the total electricity production of all PowerUtility objects in MW."""
    return sum([i.production for i in power_utilities]) / 1000 # MW

def volume_of_chemical_in_units(units, chemical):
    """Return volume of chemical that is occupied in given units [m^3]."""
    isa = isinstance
    F_vol = 0.
    for i in units:
        if isa(i, bst.BatchBioreactor):
            feed = i.ins[0]
            z_vol = feed.ivol[chemical] / feed.F_vol
            D = i.design_results
            F_vol += D['Number of reactors'] * D['Reactor volume'] * z_vol
        elif isa(i, bst.Tank):
            feed = i.ins[0]
            z_vol = feed.ivol[chemical] / feed.F_vol
            D = i.design_results
            F_vol += D['Total volume'] * z_vol
    return F_vol

def set_construction_material(units, pressure_vessel_material, tray_material,
                              tank_material, heat_exchanger_material, pump_material):
    """
    Set the construction material of all vessels, columns, trays,
    heat exchangers, and pumps
    
    Parameters
    ----------
    units : Iterable[Unit]
        Unit operations to set construction material.
    pressure_vessel_material : str
        Construction material for pressure vessels (includes distillation 
        columns and flash vessels).
    tray_material : str
        Construction material for column trays.
    tank_material : str
        Construction material for tanks.
    pump_material : str
        Construction material for pumps.
    
    """
    isa = isinstance
    HeatExchanger = bst.HX
    Distillation = bst.BinaryDistillation
    PressureVessel = bst.units.design_tools.PressureVessel
    Tank = bst.Tank
    Pump = bst.Pump
    for u in units:
        for i in (u, *u.auxiliary_units):
            if isa(i, HeatExchanger):
                i.material = heat_exchanger_material
            elif isa(i, Distillation):
                i.vessel_material = pressure_vessel_material
                i.tray_material = tray_material
            elif isa(i, PressureVessel):
                i.vessel_material = pressure_vessel_material
            elif isa(i, Tank):
                i.vessel_material = tank_material
            elif isa(i, Pump):
                i.material = pump_material

def set_construction_material_to_stainless_steel(units, kind='304'):
    """
    Set the construction material of all vessels, columns, trays,
    heat exchangers, and pumps to stainless steel.
    
    Parameters
    ----------
    kind : str, "304" or "316"
        Type of stainless steel. 
    
    """
    if kind not in ('304', '316'): 
        raise ValueError("kind must be either '304' or '316', not '%s'" %kind)
    material = 'Stainless steel ' + kind
    set_construction_material(units, material, material, 'Stainless steel',
                              'Stainless steel/stainless steel', 
                              'Stainless steel')

def set_construction_material_to_carbon_steel(units):
    """
    Set the construction material of all vessels, columns, trays,
    heat exchangers, and pumps to carbon steel or cast iron.
    
    """
    set_construction_material(units, 'Carbon steel', 'Carbon steel', 
                              'Carbon steel', 'Carbon steel/carbon steel', 
                              'Cast iron')

def get_OSBL(units):
    """
    Return a list of unit operations that would typically be classified as
    out-of-boundary-limits, just as storage tanks, cooling towers, and other
    facilities.
    
    """
    return [i for i in units
            if isinstance(i, (bst.Facility, bst.StorageTank)) and not 'surge' in i.line.lower()
            or 'storage' in i.line.lower()
            or 'wastewater' in i.line.lower()]

def heat_exchanger_operation(units):
    """
    Return a pandas DataFrame object of the capacity flow rate, 
    inlet temperature, and outlet temperature of all heat exchangers
    associated to the given units.
    
    """
    index = []
    T_in = []
    T_out = []
    C = []
    for u in units:
        for i, hu in enumerate(u.heat_utilities):
            hx = hu.heat_exchanger
            if hu.heat_exchanger:
                feed = hx.ins[0]
                index.append((u.line, u.ID, str(i)))
                T_in.append(feed.T)
                T_out.append(hx.outs[0].T)
                C.append(feed.C)
    return pd.DataFrame(
        data=np.array([C, T_in, T_out]).transpose(),
        index=pd.MultiIndex.from_tuples(index, names=['Unit', 'ID', 'index']),
        columns=['C [kJ/hr]', 'T_in [K]', 'T_out [K]'],
    )

def get_tea_results(tea, product=None, feedstock=None):
    """Return a dictionary with general TEA results."""
    TCI = tea.TCI / 1e6
    VOC = tea.VOC / 1e6
    FOC = tea.FOC / 1e6
    installed_cost = tea.installed_cost / 1e6
    material_cost = tea.material_cost / 1e6
    utility_cost = tea.utility_cost / 1e6
    sales = tea.sales / 1e6 
    dct = {
        'TCI [million USD]': round(TCI),
        'Installed equipment cost [million USD]': round(installed_cost),
        'VOC [million USD/yr]': round(VOC),
        'FOC [million USD/yr]': round(FOC),
        'Material cost [million USD/yr]': round(material_cost),
        'Utility cost [million USD/yr]': round(utility_cost),
        'Sales [million USD/yr]': round(sales),
    }
    if product:
        MPSP = tea.solve_price(product) * 907.18474
        dct['MPSP [USD/ton]'] = round(MPSP)
    if feedstock:
        MFP = tea.solve_price(feedstock) * 907.18474
        dct['MFP [USD/ton]'] = round(MFP)
    return dct
        
def default_utilities():
    """Reset utilities back to BioSTEAM's defaults."""
    bst.HeatUtility.default_agents()
    bst.PowerUtility.default_price()

def default_CEPCI():
    bst.CE = 567.5

def default(utilities=True, CEPCI=True, flowsheet=False):
    """
    Reset utilities, flowsheets, and chemical plant cost index (CEPCI) back to 
    BioSTEAM's defaults (if requested).
    
    """
    if utilities: default_utilities()
    if CEPCI: default_CEPCI()
    for i in (bst.Stream, bst.Unit, bst.System): i.ticket_numbers.clear()
    if flowsheet: 
        for i in bst.main_flowsheet.flowsheet: i.clear(False)
    bst.Chemical.chemical_cache.clear()
    
