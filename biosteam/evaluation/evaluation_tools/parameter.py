# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from chaospy import distributions as shape
import biosteam as bst

__all__ = ('load_default_parameters',
           'add_all_cost_item_params',
           'add_flow_rate_param',
           'add_all_stream_price_params',
           'add_stream_price_param',
           'add_power_utility_price_param',
           'add_basic_TEA_params',
           'triang',
           'bounded_triang',)

# %% Predefined shapes

def triang(mid, proportion=0.1, addition=0): # pragma: no cover
    lb = (1.-proportion)*mid - addition
    ub = (1.+proportion)*mid + addition
    if lb > ub: ub, lb = lb, ub
    return shape.Triangle(lb, mid, ub)

def bounded_triang(mid, lb=0, ub=1, proportion=0, addition=0.1): # pragma: no cover
    if lb > ub: ub, lb = lb, ub
    lower = (1.-proportion)*mid - addition
    upper = (1.+proportion)*mid + addition
    if lower > upper: upper, lower = lower, upper
    if lower < lb: lower = lb
    if upper > ub: upper = ub
    return shape.Triangle(lower, mid, upper)

# %% Fill model with parameters

def load_default_parameters(self, feedstock, shape=triang,
                            bounded_shape=bounded_triang,
                            operating_days=False, include_feedstock_price=True): # pragma: no cover
    """
    Load all default parameters, including coefficients of cost items,
    stream prices, electricity price, heat utility prices, feedstock
    flow rate, and number of operating days.
    
    Parameters
    ----------
    feedstock : Stream
        Main feedstock of process.
    shape : function, optional
        Should take in baseline value and return a chaospy.Dist object or
        None. Default function returns a chaospy.Triangle object with bounds
        at +- 10% of baseline value. The distribution is applied to all
        parameters except exponential factor "n" of cost items.
    bounded_shape : function, optional
        Should take in baseline value and return a chaospy.Dist object or
        None. Default function returns a chaospy.Triangle object with bounds
        at +- 0.1 of baseline value or minimum 0 and maximum 1. The 
        distribution is applied to exponential factor "n" of cost items.
    operating_days : bool, optional
        If True, include operating days
    
    """
    bounded_shape = bounded_shape or shape
    add_all_cost_item_params(self, shape, bounded_shape)
    add_all_stream_price_params(self, shape, feedstock, include_feedstock_price)
    add_power_utility_price_param(self, shape)
    add_heat_utility_price_params(self, shape)
    add_flow_rate_param(self, feedstock, shape)
    add_basic_TEA_params(self, shape, operating_days)

def add_all_cost_item_params(model, shape, exp_shape): # pragma: no cover
    system = model._system
    costunits = [i for i in system.units if i._design or i._cost]
    
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
            
def add_all_stream_price_params(model, shape, feedstock, include_feedstock_price): # pragma: no cover
    # Add stream prices as parameters
    system = model._system
    feeds = system.feeds
    products = system.products
    if not include_feedstock_price:
        feeds = feeds.copy()
        feeds.discard(feedstock)
    for s in feeds: add_stream_price_param(model, s, shape)
    for s in products: add_stream_price_param(model, s, shape)

def set_price(price): # pragma: no cover
    bst.PowerUtility.price = price
        
def add_power_utility_price_param(model, shape): # pragma: no cover
    if bst.PowerUtility.price:
        baseline = bst.PowerUtility.price
        model.parameter(set_price, element='Electricity', units='USD/kWhr',
                        distribution=shape(baseline),
                        baseline=baseline)
        
def add_heat_utility_price_params(model, shape): # pragma: no cover
    agents = (*bst.HeatUtility.cooling_agents,
                    *bst.HeatUtility.heating_agents)
    for agent in agents:
        add_agent_price_params(model, agent.ID, agent, shape)
        
class Setter: # pragma: no cover
    __slots__ = ('obj', 'attr')
    def __init__(self, obj, attr):
        self.obj = obj
        self.attr = attr
    def __call__(self, value):
        setattr(self.obj, self.attr, value)
        
def add_agent_price_params(model, name, agent, shape): # pragma: no cover
    if agent.heat_transfer_price:
        baseline = agent.heat_transfer_price
        model.parameter(Setter(agent, 'heat_transfer_price'), element=name, units='USD/kJ',
                        name='Price', distribution=shape(baseline), baseline=baseline)
    elif agent.regeneration_price:
        baseline = agent.regeneration_price
        model.parameter(Setter(agent, 'regeneration_price'),
                        element=name, units='USD/kmol',
                        name='Price', distribution=shape(baseline),
                        baseline=baseline)
        
def add_flow_rate_param(model, feed, shape): # pragma: no cover
    baseline = feed.F_mass
    model.parameter(Setter(feed, 'F_mass'), element=feed, units='kg/hr',
                    distribution=shape(baseline), kind='coupled',
                    name='Flow rate', baseline=baseline)
    
def add_basic_TEA_params(model, shape, operating_days): # pragma: no cover
    param = model.parameter
    TEA = model._system.TEA
    
    if operating_days:
        baseline = TEA.operating_days
        param(Setter(TEA, 'operating_days'), element='TEA', 
              distribution=shape(baseline),
              baseline=baseline,
              name='Operating days')
    
    baseline = TEA.income_tax
    param(Setter(TEA, 'income_tax'),
                 element='TEA',
                 distribution=shape(baseline),
                 baseline=baseline,
                 name='Income tax')
        
    if TEA.startup_months:
        baseline = TEA.startup_months
        param(Setter(TEA, 'startup_months'), element='TEA', 
              baseline=baseline,
              distribution=shape(baseline), name='startup_months')
            
    

# %% Parameter creators for cost_items

def add_stream_price_param(model, stream, shape): # pragma: no cover
    mid = stream.price
    if not mid: return
    model.parameter(Setter(stream, 'price'),
                    element=stream, units='USD/kg',
                    distribution=shape(mid),
                    baseline=mid,
                    name='price')

def _cost(model, ID, item, line, shape): # pragma: no cover
    key = 'cost'
    mid = float(item[key])
    if not mid: return None
    distribution = shape(mid)
    name = 'base cost'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    _cost_option(model, ID, item, key, line, distribution, 'USD', mid)
    
def _exp(model, ID, item, line, shape): # pragma: no cover
    key = 'n'
    mid = float(item[key])
    if not mid: return None
    distribution = shape(mid)
    name = 'exponent'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    _cost_option(model, ID, item, key, line, distribution, None, mid)

class CostItemSetter: # pragma: no cover
    __slots__ = ('item', 'key', 'size')
    def __init__(self, item, key, size=1):
        self.item = item
        self.key = key
        self.size = size
    def __call__(self, value):
        self.item[self.key] = value*self.size

def _kW(model, ID, item, line, shape): # pragma: no cover
    key = 'kW'
    mid = float(item[key])
    if not mid: return None
    units = item.units
    size = item.S
    mid = mid/size
    distribution = shape(mid)
    units = 'kW / ' + units
    name = 'electricity rate'
    if ID!=line: ID = ID + ' ' + name
    else: ID = name
    model.parameter(CostItemSetter(item, key, size), element=line,
                    baseline=mid,
                    units=units, distribution=distribution, name=ID)

def _cost_option(model, ID, item, key, line, distribution, units, mid): # pragma: no cover
    model.parameter(CostItemSetter(item, key), element=line, units=units, 
                    baseline=mid,
                    distribution=distribution, name=ID)
 