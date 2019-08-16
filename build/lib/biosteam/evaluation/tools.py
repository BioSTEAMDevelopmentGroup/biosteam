# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:21:31 2019

@author: yoelr
"""
from chaospy import distributions as shape
import biosteam as bst
import pandas as pd

__all__ = ('add_all_parameters',
           'add_all_cost_item_params',
           'add_flow_rate_param',
           'add_all_stream_price_params',
           'add_stream_price_param',
           'add_power_utility_price_param',
           'add_basic_TEA_params',
           'triang',
           'bounded_triang',)

# %% Predefined shapes

def triang(mid, proportion=0.1, addition=0):
    return shape.Triangle((1.-proportion)*mid - addition,
                          mid,
                          (1.+proportion)*mid + addition)

def bounded_triang(mid, lb=0, ub=1, proportion=0, addition=0.1):
    lower = (1.-proportion)*mid - addition
    upper = (1.+proportion)*mid + addition
    if lower < lb: lower = lb
    if upper > ub: upper = ub
    return shape.Triangle(lower, mid, upper)

# %% Fill model with parameters

def add_all_parameters(model, feedstock, shape=triang,
                       bounded_shape=bounded_triang):
    bounded_shape = bounded_shape or shape
    add_all_cost_item_params(model, shape, bounded_shape)
    add_all_stream_price_params(model, shape)
    add_power_utility_price_param(model, shape)
    add_heat_utility_price_params(model, shape)
    add_flow_rate_param(model, feedstock, shape)
    add_basic_TEA_params(model, shape)

def add_all_cost_item_params(model, shape, exp_shape):
    system = model._system
    costunits = system._costunits
    
    # Get all cost options (without repetition)
    all_unit_lines = []
    all_cost_items = []
    for i in costunits:
        if hasattr(i, 'cost_items'):
            if any([i.cost_items is j for j in all_cost_items]): continue
            all_unit_lines.append(i.line)    
            all_cost_items.append(i.cost_items)
    
    # Add cost exponents, base cost, and electricity cost as parameters
    # (for each column in each dataframe of the list)
    for cost_items, line in zip(all_cost_items, all_unit_lines):
        for ID, item in cost_items.items():
            _cost(model, ID, item, line, shape)
            _exp(model, ID, item, line, exp_shape)
            _kW(model, ID, item, line, shape)
            
def add_all_stream_price_params(model, shape):
    # Add stream prices as parameters
    system = model._system
    for s in system.feeds: add_stream_price_param(model, s, shape)
    for s in system.products: add_stream_price_param(model, s, shape)
    
def add_power_utility_price_param(model, shape):
    if bst.PowerUtility.price:
        @model.parameter(element='Electricity', units='USD/kW',
                         distribution=shape(bst.PowerUtility.price))
        def set_price(price):
            bst.PowerUtility.price = price 

def add_heat_utility_price_params(model, shape):
    named_agents = (*bst.HeatUtility.cooling_agents.items(),
                    *bst.HeatUtility.heating_agents.items())
    for name, agent in named_agents:
        add_agent_price_params(model, name, agent, shape)
        
def add_agent_price_params(model, name, agent, shape):
    if agent.price_kJ:
        @model.parameter(element=name, units='USD/kJ',
                         name='Price per unit energy',
                         distribution=shape(agent.price_kJ))
        def set_price(price_kJ):
            agent.price_kJ = price_kJ 
    elif agent.price_kmol:
        @model.parameter(element=name, units='USD/kmol',
                         name='Price per unit material',
                         distribution=shape(agent.price_kmol))
        def set_price(price_kmol):
            agent.price_kmol = price_kmol 
        
def add_flow_rate_param(model, feed, shape):
    @model.parameter(element=feed, units='kg/hr',
                     distribution=shape(feed.massnet),
                     kind='coupled')
    def set_flow_rate(flow_rate):
        feed.mol[:] *= flow_rate/feed.massnet
    
def add_basic_TEA_params(model, shape):
    param = model.parameter
    TEA = model._system.TEA
    
    @param(element='TEA', distribution=shape(TEA.operating_days))
    def set_operating_days(operating_days):
        TEA.operating_days = operating_days
        
    @param(element='TEA', distribution=shape(TEA.income_tax))
    def set_income_tax(income_tax):
        TEA.income_tax = income_tax
        
    @param(element='TEA', distribution=shape(TEA.startup_months))
    def set_startup_months(startup_months):
        TEA.startup_months = startup_months
        
    

# %% Parameter creators for cost_items

def add_stream_price_param(model, stream, shape):
    mid = stream.price
    if not mid: return
    @model.parameter(element=stream, units='USD/kg', distribution=shape(mid))
    def set_price(price): stream.price = price

def _cost(model, ID, item, line, shape):
    key = 'cost'
    mid = float(item[key])
    distribution = shape(mid)
    name = 'base cost'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    _cost_option(model, ID, item, key, line, distribution, 'USD')
    
def _exp(model, ID, item, line, shape):
    key = 'n'
    mid = float(item[key])
    distribution = shape(mid)
    name = 'exponent'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    _cost_option(model, ID, item, key, line, distribution, None)

def _kW(model, ID, item, line, shape):
    key = 'kW'
    mid = float(item[key])
    if not mid: return None
    distribution = shape(mid)
    basis = item.basis
    units = item.units
    size = item.S
    size = '' if size is 1 else f'{size} '
    units = f'kW per {size}{units} {basis.casefold()}'
    name = 'electricity rate'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    _cost_option(model, ID, item, key, line, distribution, units)

def _cost_option(model, ID, item, key, line, distribution, units):
    @model.parameter(element=line, units=units, distribution=distribution, name=ID)
    def set_cost_option(value):
        item[key] = value
    set_cost_option
 