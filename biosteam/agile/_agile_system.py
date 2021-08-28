# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains the abstract AgileSystem class that serves to create
objects that can retrive general results from multiple system scenarios
in such a way that it represents an agile production process.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.agile.AgileSystem

"""

from . import AgileScenario
from .._system import ScenarioCosts
from inspect import signature
from thermosteam.utils import repr_kwargs

__all__ = ('AgileSystem',)

class OperationMode:
    __slots__ = ('operating_hours', 'data')
    def __init__(self, operating_hours, **data):
        self.operating_hours = operating_hours
        self.data = data
    
    def __getattr__(self, name):
        if name in self.data:
            return self.data[name]
        elif name in self.agile_system.operation_parameters:
            raise AttributeError(f"operation parameter '{name}' is defined, but has no value")
        else:
            raise AttributeError(f"'{name}' is not a defined operation parameter")
    
    def __setattr__(self, name, value):
        if name in self.agile_system.operation_parameters:
            self.data[name] = value
        elif name == 'operating_hours':
            object.__setattr__(self, 'operating_hours', value)
        elif name == 'data':
            object.__setattr__(self, 'data', value)
        else:
            raise AttributeError(f"'{name}' is not a defined operation parameter")

    def __repr__(self):
        return f"{type(self).__name__}(operating_hours={self.operating_hours}{repr_kwargs(self.data)})"
    
        
class AgileSystem(AgileScenario):
    """
    Class for creating objects which may serve to retrive
    general results from multiple system scenarios in such a way that it 
    represents an agile production process. When simulated, an AgileSystem 
    creates scenarios from simulated system operation modes and compile them to 
    retrieve results later.
    
    Parameters
    ----------
    system : System
        System object that is simulated in each operation mode.
    operation_modes : list[OperationMode], optional
        Defines each mode of operation with time steps and parameter values
    operation_parameters : dict[str: function], optional
        Defines all parameters available for all operation modes.
    tea : AgileTea, optional
        If given, the AgileTEA object will compile results after each simulation.
    
    
    """
    
    __slots__ = ('system', 'operation_modes', 'operation_parameters',
                 'unit_capital_costs', 'utility_cost',
                 'flow_rates', 'feeds', 'products',
                 'tea', '_OperationMode')
    
    def __init__(self, system, operation_modes=None, operation_parameters=None, tea=None):
        self.system = system
        self.operation_modes = [] if operation_modes is None else operation_modes 
        self.operation_parameters = {} if operation_parameters  is None else operation_parameters 
        self.tea = tea
        self._OperationMode = type('OperationMode', (OperationMode,), {'agile_system': self})
        
    def _downstream_system(self, unit):
        return self

    def operation_mode(self, operating_hours, **data):
        """
        Define and register an operation mode.
        
        Parameters
        ----------    
        operating_hours : function
            Length of operation in hours.
        **data : str
            Name and value-pairs of operation parameters.
        
        """
        om = self._OperationMode(operating_hours, **data)
        self.operation_modes.append(om)
        return om

    def operation_parameter(self, setter, name=None):
        """
        Define and register operation parameter.
        
        Parameters
        ----------    
        setter : function
                 Should set parameter in the element.
        name : str
               Name of parameter. If None, default to argument name of setter.
        
        """
        if not setter: return lambda setter: self.operation_parameter(setter, name)
        if not name: name, *_ = signature(setter).parameters.keys()
        self.operation_parameters[name] = setter
        return setter

    @property
    def units(self):
        return self.system.units

    @property
    def empty_recycles(self):
        return self.system.empty_recycles
    
    @property    
    def reset_cache(self):
        return self.system.reset_cache

    @property
    def operating_hours(self):
        return sum([i.operating_hours for i in self.operation_modes])

    def create_scenario(self, system):
        return system.get_scenario_costs()
            
    def compile_scenarios(self, scenarios):
        units = set(sum([list(i.unit_capital_costs) for i in scenarios], []))
        unit_scenarios = {i: [] for i in units}
        for scenario in scenarios:
            unit_capital_costs = scenario.unit_capital_costs
            for i, j in unit_capital_costs.items(): unit_scenarios[i].append(j)
        self.unit_capital_costs = {i: i.get_agile_capital_costs(j) for i, j in unit_scenarios.items()}
        self.utility_cost = sum([i.utility_cost for i in scenarios])
        self.flow_rates = flow_rates = {}
        self.feeds = list(set(sum([i.feeds for i in scenarios], [])))
        self.products = list(set(sum([i.products for i in scenarios], [])))
        for scenario in scenarios:
            for stream, F_mass in scenario.flow_rates.items():
                if stream in flow_rates: flow_rates[stream] += F_mass
                else: flow_rates[stream] = F_mass
    
    def simulate(self):
        scenarios = []
        system = self.system
        original_operating_hours = system.operating_hours
        for operation_mode in self.operation_modes:
            operating_hours = operation_mode.operating_hours
            if operating_hours == 0.: continue
            operation_parameters = self.operation_parameters
            for name, value in operation_mode.data.items():    
                operation_parameters[name](value)
            system.operating_hours = operating_hours
            system.simulate()
            scenario = self.create_scenario(system)
            scenarios.append(scenario)
        system.operating_hours = original_operating_hours
        self.compile_scenarios(scenarios)
        if self.tea:
            scenario = self.tea.create_scenario(self)
            self.tea.compile_scenarios([scenario])
        
    def get_scenario_costs(self):
        return ScenarioCosts(
            self.unit_capital_costs,
            self.flow_rates,
            self.utility_cost,
            self.feeds, self.products,
            self.operating_hours,
        )
    
            
