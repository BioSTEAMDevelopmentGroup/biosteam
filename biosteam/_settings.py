# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module extends thermosteam's ProcessSettings object.
"""
from thermosteam import settings
from thermosteam.units_of_measure import UnitsOfMeasure
from typing import Callable
import biosteam as bst
import numpy as np

# %% Excecutables for dynamic programing of new stream utilities

get_unit_utility_flow_executable = '''
def {flowmethname}(self):
    """Return the {docname} flow rate [kg/hr]."""
    flow = 0.
    if {name} in self._inlet_utility_indices:
        flow += self._ins[self._inlet_utility_indices[{name}]].F_mass
    if {name} in self._outlet_utility_indices:
        flow += self._outs[self._outlet_utility_indices[{name}]].F_mass
    return flow

bst.Unit.{flowmethname} = {flowmethname}
'''

get_unit_utility_cost_executable = '''
def {costmethname}(self):
    """Return the {docname} cost [USD/hr]."""
    return bst.stream_prices[{name}] * self.{flowmethname}()

bst.Unit.{costmethname} = {costmethname}
'''

get_system_utility_flow_executable = '''
def {flowmethname}(self):
    """Return the {docname} flow rate [kg/yr]."""
    return sum([i.get_utility_flow({name}) for i in self.cost_units]) * self.operating_hours

bst.System.{flowmethname} = {flowmethname}
'''

get_system_utility_cost_executable = '''
def {costmethname}(self):
    """Return the {docname} cost [USD/yr]."""
    return bst.stream_prices[name] * self.{flowmethname}()

bst.System.{costmethname} = {costmethname}
'''

# %% New properties

@property
def CEPCI(self) -> float:
    """Chemical engineering plant cost index (defaults to 567.5 at 2017)."""
    return bst.CE
@CEPCI.setter
def CEPCI(self, CEPCI):
    bst.CE = CEPCI

@property
def utility_characterization_factors(self) ->  dict[tuple[str, str], tuple[float, UnitsOfMeasure]]:
    """Utility characterization factor data (value and units) by agent ID 
    and impact key."""
    return bst.HeatUtility.characterization_factors
@utility_characterization_factors.setter
def utility_characterization_factors(self, utility_characterization_factors):
    bst.HeatUtility.characterization_factors = utility_characterization_factors

@property
def cooling_agents(self) -> list[bst.UtilityAgent]:
    """All cooling utilities available."""
    return bst.HeatUtility.cooling_agents
@cooling_agents.setter
def cooling_agents(self, cooling_agents):
    bst.HeatUtility.cooling_agents = cooling_agents
    
@property
def heating_agents(self) -> list[bst.UtilityAgent]:
    """All heating utilities available."""
    return bst.HeatUtility.heating_agents
@heating_agents.setter
def heating_agents(self, heating_agents):
    bst.HeatUtility.heating_agents = heating_agents
    
@property
def stream_prices(self) -> dict[str, float]:
    """Price of stream utilities/fees/credits [USD/kg] which are defined as 
    inlets and outlets to unit operations."""
    return bst.stream_prices
@stream_prices.setter
def stream_prices(self, stream_prices):
    bst.stream_prices = stream_prices

@property
def impact_indicators(self) -> dict[str, str]:
    """User-defined impact indicators and their units of measure."""
    return bst.impact_indicators
@impact_indicators.setter
def impact_indicators(self, impact_indicators):
    bst.impact_indicators = impact_indicators

@property
def electricity_price(self) -> float:
    """Electricity price [USD/kWhr]"""
    return bst.PowerUtility.price
@electricity_price.setter
def electricity_price(self, electricity_price):
    """Electricity price [USD/kWhr]"""
    bst.PowerUtility.price = electricity_price

@property
def skip_simulation_of_units_with_empty_inlets(self):
    """Whether to assume unit does not exist when all inlet streams are empty. 
    If inlets are empty and this flag is True, detailed mass and energy 
    balance, design, and costing algorithms are skipped and all outlet streams 
    are emptied."""
    return bst.Unit._skip_simulation_when_inlets_are_empty
@skip_simulation_of_units_with_empty_inlets.setter
def skip_simulation_of_units_with_empty_inlets(self, skip):
    bst.Unit._skip_simulation_when_inlets_are_empty = skip

@property
def allocation_properties(self):
    """Defined allocation property and basis pairs for LCA."""
    return bst.allocation_properties

def register_utility(self, name: str, price: float):
    """Register new stream utility/credit/fee in BioSTEAM given the name and the price 
    [USD/kg]."""
    if name not in bst.stream_prices:
        docname = name.lower()
        methname = docname.replace(' ', '_').replace('-', '_')
        flowmethname = f"get_{methname}_flow"
        costmethname = f"get_{methname}_cost"
        repname = repr(name)
        globs = {'bst': bst}
        flow_kwargs = dict(
            flowmethname=flowmethname,
            docname=docname,
            name=repname,
        )
        cost_kwargs = dict(
            costmethname=costmethname,
            flowmethname=flowmethname,
            docname=docname,
            name=repname,
        )
        
        # Unit
        exec(get_unit_utility_flow_executable.format(**flow_kwargs), globs)
        exec(get_unit_utility_cost_executable.format(**cost_kwargs), globs)
        
        # System
        exec(get_system_utility_flow_executable.format(**flow_kwargs), globs)
        exec(get_system_utility_cost_executable.format(**cost_kwargs), globs)
    bst.stream_prices[name] = price

def define_allocation_property(
        self, name: str, basis: float, 
        stream: Callable=None, 
        power_utility: Callable=None,
        heat_utility: Callable=None
    ):
    """Define a new allocation property by property getters."""
    allocation_name = name + '-allocation'
    units = basis + '/hr'
    if stream is not None:
        bst.Stream.define_property(
            allocation_name, units, stream,
        )
    if power_utility is not None:
        bst.PowerUtility.define_property(
            allocation_name, units, power_utility,
        )
    if heat_utility is not None:
        bst.HeatUtility.define_property(
            allocation_name, units, heat_utility,
        )
    bst.allocation_properties[name] = basis

Settings = settings.__class__
Settings.CEPCI = CEPCI
Settings.utility_characterization_factors = utility_characterization_factors
Settings.cooling_agents = cooling_agents
Settings.heating_agents = heating_agents
Settings.impact_indicators = impact_indicators
Settings.stream_prices = stream_prices
Settings.electricity_price = electricity_price
Settings.skip_simulation_of_units_with_empty_inlets = skip_simulation_of_units_with_empty_inlets
Settings.register_fee = Settings.register_credit = Settings.register_utility = register_utility
Settings.allocation_properties = allocation_properties
Settings.define_allocation_property = define_allocation_property

# %% Register stream utilities

settings.register_utility('Fuel', 0.218) # Natural gas
settings.register_utility('Natural gas', 0.218) # Natural gas
settings.register_utility('Ash disposal', -0.0318)
settings.register_utility('Reverse osmosis water', 5.6e-4)
settings.register_utility('Process water', 2.7e-4)

# %% Register predefined allocation methods

settings.define_allocation_property(
    'energy', 'kJ', 
    stream=lambda self: max(self.LHV, 0),
    power_utility=lambda self: max(-self.rate * 3600, 0.),
    heat_utility=lambda self: max(self.duty, 0),
)
settings.define_allocation_property(
    'revenue', 'USD', 
    stream=lambda self: self.cost if self.price > 0. else 0.,
    power_utility=lambda self: max(-self.cost, 0.),
    heat_utility=lambda self: max(self.cost, 0),
)
settings.define_allocation_property(
    'mass', 'kg', stream=lambda self: np.dot(self.chemicals.MW, self.mol)
)
