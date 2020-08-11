# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import numpy as np
import biosteam as bst
from .._heat_utility import HeatUtility

__all__ = ('heat_exchanger_utilities_from_units',
           'stream_mass_balance',
           'units_with_costs',
           'get_utility_flow',
           'get_utility_duty',
           'get_power_utilities',
           'get_heat_utilities',
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
           'default_utilities',
)
    
def units_with_costs(units):
    """Filter units that have a cost or design method."""
    return [i for i in units if i._cost or i._design]

def heat_exchanger_utilities_from_units(units):
    """Return a list of heat utilities from all heat exchangers,
    including the condensers and boilers of distillation columns and
    flash vessel heat exchangers."""
    heat_utilities = sum([i.heat_utilities for i in units], ())
    return [i for i in heat_utilities if i.heat_exchanger]

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

def get_utility_flow(heat_utilities, agent):
    """Return the total utility duty of heat utilities for given agent in GJ/hr"""
    if isinstance(agent, str):
        agent = HeatUtility.get_agent(agent)
    return sum([i.flow * i.agent.MW for i in heat_utilities if i.duty > 0. and i.agent is agent]) / 1e3

def get_utility_duty(heat_utilities, agent):
    """Return the total utility duty of heat utilities for given agent in GJ/hr"""
    if isinstance(agent, str):
        agent = HeatUtility.get_agent(agent)
    return sum([i.duty for i in heat_utilities if i.duty > 0. and i.agent is agent]) / 1e6 

def get_power_utilities(units):
    """Return a list of all PowerUtility objects."""
    return [i.power_utility for i in units if i.power_utility]

def get_heat_utilities(units):
    """Return a list of all HeatUtility objects."""
    return sum([i.heat_utilities for i in units], ())

def get_installed_cost(units):
    """Return the total installed equipment cost of all units in million USD."""
    return sum([i.installed_cost for i in units]) / 1e6 # millions USD

def get_cooling_duty(heat_utilities):
    """Return the total cooling duty of all heat utilities in GJ/hr."""
    cooling_duty = sum([i.duty for i in heat_utilities if i.duty < 0 and i.flow > 0]) / 1e6 # GJ/hr
    return abs(cooling_duty)
               
def get_heating_duty(heat_utilities):
    """Return the total heating duty of all heat utilities in GJ/hr."""
    return sum([i.duty for i in heat_utilities if i.duty > 0 and i.flow > 0]) / 1e6 # GJ/hr
          
def get_electricity_consumption(power_utilities):
    """Return the total electricity consumption of all PowerUtility objects in MW."""
    return sum([i.consumption for i in power_utilities]) / 1000 # MW

def get_electricity_production(power_utilities):
    """Return the total electricity production of all PowerUtility objects in MW."""
    return sum([i.production for i in power_utilities]) / 1000 # MW

def volume_of_chemical_in_units(units, chemical):
    """Return volume of chemical that occupied in given units [m^3]."""
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
    for i in units:
        if isa(i, HeatExchanger):
            i.material = heat_exchanger_material
        elif isa(i, Distillation):
            i.boiler.material = i.condenser.material = heat_exchanger_material
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

def stream_mass_balance(chemical_IDs, variable_inlets, constant_inlets=(),
                        constant_outlets=(), is_exact=True, balance='flow'):
    """
    Solve stream mass balance by iteration.
    
    Parameters
    ----------
    chemical_IDs : tuple[str]
        Chemicals that will be used to solve mass balance linear equations.
        The number of chemicals must be same as the number of input streams varied.
    variable_inlets : Iterable[Stream]
        Inlet streams that can vary in net flow rate to accomodate for the
        mass balance.
    constant_inlets: Iterable[Stream], optional
        Inlet streams that cannot vary in flow rates.
    constant_outlets: Iterable[Stream], optional
        Outlet streams that cannot vary in flow rates.
    is_exact=True : bool, optional
        True if exact flow rate solution is required for the specified IDs.
    kind='flow' : {'flow', 'composition'}, optional
          * 'flow': Satisfy output flow rates
          * 'composition': Satisfy net output molar composition
    
    """
    # SOLVING BY ITERATION TAKES 15 LOOPS FOR 2 STREAMS
    # SOLVING BY LEAST-SQUARES TAKES 40 LOOPS
    solver = np.linalg.solve if is_exact else np.linalg.lstsq

    # Set up constant and variable streams
    if not variable_inlets:
        raise ValueError('variable_inlets must contain at least one stream')
    index = variable_inlets[0].chemicals.get_index(chemical_IDs)
    mol_out = sum([s.mol for s in constant_outlets])

    if balance == 'flow':
        # Perform the following calculation: Ax = b = f - g
        # Where:
        #    A = flow rate array
        #    x = factors
        #    b = target flow rates
        #    f = output flow rates
        #    g = constant inlet flow rates

        # Solve linear equations for mass balance
        A = np.array([s.mol for s in variable_inlets]).transpose()[index, :]
        f = mol_out[index]
        g = sum([s.mol[index] for s in constant_inlets])
        b = f - g
        x = solver(A, b)

        # Set flow rates for input streams
        for factor, s in zip(x, variable_inlets):
            s.mol[:] = s.mol * factor

    elif balance == 'composition':
        # Perform the following calculation:
        # Ax = b
        #    = sum( A_ * x_guess + g_ )f - g
        #    = A_ * x_guess * f - O
        # O  = sum(g_)*f - g
        # Where:
        # A_ is flow array for all species
        # g_ is constant flows for all species
        # Same variable definitions as in 'flow'

        # Set all variables
        A_ = np.array([s.mol for s in variable_inlets]).transpose()
        A = np.array([s.mol for s in variable_inlets]).transpose()[index, :]
        F_mol_out = mol_out.sum()
        z_mol_out = mol_out / F_mol_out if F_mol_out else mol_out
        f = z_mol_out[index]
        g_ = sum([s.mol for s in constant_inlets])
        g = g_[index]
        O = sum(g_) * f - g

        # Solve by iteration
        x_guess = np.ones_like(index)
        not_converged = True
        while not_converged:
            # Solve linear equations for mass balance
            b = (A_ * x_guess).sum()*f + O
            x_new = solver(A, b)
            not_converged = sum(((x_new - x_guess)/x_new)**2) > 0.0001
            x_guess = x_new

        # Set flow rates for input streams
        for factor, s in zip(x_new, variable_inlets):
            s.mol = s.mol * factor
    
    else:
        raise ValueError( "balance must be one of the following: 'flow', 'composition'")
        
def default_utilities():
    """Reset utilities back to BioSTEAM's defaults."""
    bst.HeatUtility.default_agents()
    bst.PowerUtility.default_price()