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

__all__ = ('AgileSystem',)

class AgileSystem(AgileScenario):
    """
    Abstract class for creating objects which may serve to retrive
    general results from multiple system scenarios in such a way that it 
    represents an agile production process. When simulated, an AgileSystem 
    creates scenarios from simulated system samples and compile them to 
    retrieve results later.
    
    Abstract Methods
    ----------------
    set_parameters(samples)
    
    """
    
    __slots__ = ('system', 'samples', 'time_steps',
                 'unit_capital_costs', 'utility_cost',
                 'flow_rates', 'feeds', 'products',
                 'tea')
    
    def __init__(self, system, samples, time_steps, tea=None):
        self.system = system
        self.samples = samples
        self.time_steps = time_steps
        self.tea = tea
    
    def __init_subclass__(cls):
        if not hasattr(cls, 'set_parameters'):
            raise NotImplementedError("missing method 'set_parameters'")

    def _downstream_system(self, unit):
        return self

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
        return sum(self.time_steps)

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
        for operating_hours, sample in zip(self.time_steps, self.samples):
            self.set_parameters(sample)
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
            sum(self.time_steps),
        )
    
            
