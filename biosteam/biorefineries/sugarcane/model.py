# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:21:31 2019

@author: yoelr
"""
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import triang
import biosteam.biorefineries.sugarcane as sc

__all__ = ('sugarcane_model',)

tea = sc.sugarcane_sys.TEA
ethanol = sc.system.ethanol
products = (ethanol,)
get_prodcost = lambda: float(tea.production_cost(products))
get_FCI = lambda: tea._FCI_cached
get_prod = lambda: ethanol.massnet * tea._annual_factor

BT = sc.system.BT
sc_sys = sc.sugarcane_sys
def get_steam():
    return sum([i.flow for i in BT.steam_utilities])*18.01528*tea._annual_factor/1000

power_utils = ([i._power_utility for i in sc_sys.units if (i._power_utility and i is not BT)])
excess_electricity = [0]
def get_consumed_electricity():
    factor = tea._annual_factor/1000
    electricity_generated = -BT._power_utility.rate * factor
    consumed_electricity = sum([i.rate for i in power_utils]) * factor
    excess_electricity[0] = electricity_generated - consumed_electricity
    return consumed_electricity
get_excess_electricity = lambda: excess_electricity[0]

metrics = (Metric('Internal rate of return', sc.sugarcane_tea.solve_IRR, '%'),
           Metric('Ethanol production cost', get_prodcost, 'USD/yr'),
           Metric('Fixed capital investment', get_FCI, 'USD'),
           Metric('Ethanol production', get_prod, 'kg/hr'),
           Metric('Steam', get_steam, 'MT/yr'),
           Metric('Consumed electricity', get_consumed_electricity, 'MWhr/yr'),
           Metric('Excess electricity', get_excess_electricity, 'MWhr/yr'))

sugarcane_model = Model(sc.sugarcane_sys, metrics, skip=False)
sugarcane_model.load_default_parameters(sc.system.Sugar_cane)
param = sugarcane_model.parameter

# Fermentation efficiency
fermentation = sc.system.P24
@param(element=fermentation, distribution=triang(fermentation.efficiency),
       kind='coupled')
def set_fermentation_efficiency(efficiency):
    fermentation.efficiency= efficiency
    
# Boiler efficiency
BT = sc.system.BT
@param(element=BT, distribution=triang(BT.boiler_efficiency))
def set_boiler_efficiency(boiler_efficiency):
    BT.boiler_efficiency = boiler_efficiency

# Turbogenerator efficiency
@param(element=BT, distribution=triang(BT.turbogenerator_efficiency))
def set_turbogenerator_efficiency(turbo_generator_efficiency):
    BT.turbo_generator_efficiency = turbo_generator_efficiency
    
# RVF separation
rvf = sc.system.P14
@param(element=rvf, distribution=triang(rvf.split['Lignin']),
        kind='coupled')
def set_rvf_solids_retention(solids_retention):
    rvf.split['Lignin', 'CaO', 'Ash', 'Cellulose', 'Hemicellulose'] = solids_retention









