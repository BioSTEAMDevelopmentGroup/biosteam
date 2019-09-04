# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:43:56 2019

@author: yoelr
"""
from biosteam.evaluation import evaluation_tools as tools
from biosteam.evaluation import Model, Metric
from biosteam.biorefineries.cornstover.system import \
    cornstover_sys, ethanol_tea, cornstover_tea, \
    ethanol, feed, F300, ethanol_density_kggal

get_MESP = lambda: cornstover_tea.solve_price(ethanol, ethanol_tea) * ethanol_density_kggal

get_FCI = lambda: sum([i._FCI_cached for i in cornstover_tea.TEAs])

get_coproduct_credit = lambda: sum([i._utility_cost_cached for i in cornstover_tea.TEAs])

get_ethanol_sales = lambda: ethanol.cost * ethanol_tea._annual_factor

metrics =(Metric('Minimum ethanol selling price', 'USD/gal', get_MESP),
          Metric('Fixed capital investment', 'USD', get_FCI),
          Metric('Co-product credit', 'USD/yr', get_coproduct_credit),
          Metric('Ethanol sales', 'USD/yr', get_ethanol_sales))

cornstover_model = Model(cornstover_sys, metrics)
cornstover_model.load_default_parameters(feed)
param = cornstover_model.parameter

# Add saccharification as a parameter
saccharification_reaction = F300.saccharification[2]
X = tools.bounded_triang(saccharification_reaction.X, addition=0.05)
@param(element=F300, kind='coupled', distribution=X)
def set_saccharification_conversion(saccharification_conversion):
    saccharification_reaction.X = saccharification_conversion

# Add ethanol conversion as a parameter
ethanol_reaction = F300.cofermentation[0]
X = tools.bounded_triang(ethanol_reaction.X, addition=0.05)
@param(element=F300, kind='coupled', distribution=X)
def set_ethanol_conversion(ethanol_conversion):
    ethanol_reaction.X = ethanol_conversion
    
# Add saccharification time as a parameter
X = tools.triang(F300.tau_saccharification)
@param(element=F300, kind='isolated', distribution=X)
def set_saccharification_time(saccharification_time):
    F300.tau_saccharification= saccharification_time
    
# Add fermentation time as a parameter
X = tools.triang(F300.tau_cofermentation)
@param(element=F300, kind='isolated',  distribution=X)
def set_fermentation_time(fermentation_time):
    F300.tau_cofermentation = fermentation_time

