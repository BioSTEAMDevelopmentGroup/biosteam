# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:43:56 2019

@author: yoelr
"""
from biosteam.evaluation import evaluation_tools as tools
from biosteam.evaluation import Model, Metric
from biosteam.biorefineries.cornstover.system import \
    cornstover_sys, ethanol_tea, cornstover_tea, \
    ethanol, cornstover, R301, ethanol_density_kggal, \
    areas, BT, Area700

get_MESP = lambda: cornstover_tea.solve_price(ethanol, ethanol_tea) * ethanol_density_kggal
get_FCI = lambda: sum([i._FCI_cached for i in cornstover_tea.TEAs])
get_coproduct_credit = lambda: sum([i._utility_cost_cached for i in cornstover_tea.TEAs])
get_ethanol_production = lambda: ethanol.massnet
get_steam_demand = lambda: BT.steam_demand.massnet
pws = [i._power_utility for i in cornstover_sys.units
       if i._power_utility and i._power_utility.rate > 0]
get_excess_electricity = lambda: BT._Design['Work'] - sum([i.rate for i in pws])

metrics =[Metric('Minimum ethanol selling price', get_MESP, 'USD/gal'),
          Metric('Fixed capital investment', get_FCI, 'USD'),
          Metric('Co-product credit', get_coproduct_credit, 'USD/yr'),
          Metric('Ethanol production', get_ethanol_production, 'kg/hr'),
          Metric('Steam demand', get_steam_demand, 'kg/hr'),
          Metric('Excess electricity', get_excess_electricity, 'kW')]

def electricity_rate_function(tea):
    power_utilities = [i._power_utility for i in tea.units if i._has_power_utility]
    if tea is Area700:
        boiler_item = BT.cost_items['Boiler']
        Design = BT._Design
        return lambda: boiler_item.kW/boiler_item.S * Design['Flow rate']/1e3
    return lambda: sum([i.rate for i in power_utilities])/1e3

def cooling_duty_function(tea):
    heat_utilities = sum([i._heat_utilities for i in tea.units if i._N_heat_utilities], [])
    cooling_utilities = [i for i in heat_utilities if i.duty < 0]
    return lambda: sum([i.duty for i in cooling_utilities])

def installation_cost_function(tea):
    return lambda: tea.installation_cost

for i, tea in enumerate(areas, 1):
    Area = f'Area {i}00'
    metrics.extend(
        (Metric('Electricity', electricity_rate_function(tea), 'MW', Area),
         Metric('Cooling duty', cooling_duty_function(tea), 'MMkcal/hr', Area),
         Metric('Installation cost', installation_cost_function(tea), '10^6 USD', Area)))

cornstover_model = Model(cornstover_sys, metrics)
cornstover_model.load_default_parameters(cornstover, operating_days=False)
param = cornstover_model.parameter

# Add saccharification as a parameter
saccharification_reaction = R301.saccharification[2]
X = tools.bounded_triang(saccharification_reaction.X, addition=0.05)
@param(element=R301, kind='coupled', distribution=X)
def set_saccharification_conversion(saccharification_conversion):
    saccharification_reaction.X = saccharification_conversion

# Add ethanol conversion as a parameter
ethanol_reaction = R301.cofermentation[0]
X = tools.bounded_triang(ethanol_reaction.X, addition=0.05)
@param(element=R301, kind='coupled', distribution=X)
def set_ethanol_conversion(ethanol_conversion):
    ethanol_reaction.X = ethanol_conversion
    
# Add saccharification time as a parameter
X = tools.triang(R301.tau_saccharification)
@param(element=R301, kind='isolated', distribution=X)
def set_saccharification_time(saccharification_time):
    R301.tau_saccharification= saccharification_time
    
# Add fermentation time as a parameter
X = tools.triang(R301.tau_cofermentation)
@param(element=R301, kind='isolated',  distribution=X)
def set_fermentation_time(fermentation_time):
    R301.tau_cofermentation = fermentation_time

