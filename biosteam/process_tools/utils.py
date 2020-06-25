# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from .._heat_utility import HeatUtility

__all__ = ('get_utility_flow',
           'get_utility_duty',
           'get_power_utilities',
           'get_heat_utilities',
           'get_installed_cost',
           'get_cooling_duty',
           'get_heating_duty',
           'get_electricity_consumption',
           'get_electricity_production',
)

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