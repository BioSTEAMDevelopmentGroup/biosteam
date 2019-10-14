# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:21:31 2019

@author: yoelr
"""
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import triang
import biosteam.biorefineries.lipidcane as lc
import numpy as np

__all__ = ('lipidcane_model', 'lipidcane_model_with_lipidfraction_parameter')

tea = lc.lipidcane_tea
ethanol = lc.system.ethanol.ethanol
biodiesel = lc.system.biodiesel.biodiesel
lipid_cane = lc.system.pretreatment.lipid_cane

etoh_prodcost = [0]
products = (biodiesel, ethanol)
def get_biodiesel_prodcost():
    bd, etoh_prodcost[0] = tea.production_cost(products)
    return bd
get_etoh_prodcost = lambda: etoh_prodcost[0]
get_FCI = lambda: tea._FCI_cached

etoh_prod = [0]
def get_biodiesel_prod():
    bd, etoh_prod[0] = np.array([biodiesel.massnet, ethanol.massnet]) * tea._annual_factor
    return bd
get_etoh_prod = lambda: etoh_prod[0]

BT = lc.system.biorefinery.BT
lc_sys = lc.lipidcane_sys
def get_steam():
    return sum([i.flow for i in BT.steam_utilities])*18.01528*tea._annual_factor/1000

power_utils = ([i._power_utility for i in lc_sys.units if (i._power_utility and i is not BT)])
excess_electricity = [0]
def get_consumed_electricity():
    factor =  tea._annual_factor/1000
    electricity_generated = -BT._power_utility.rate * factor
    consumed_electricity = sum([i.rate for i in power_utils]) * factor
    excess_electricity[0] = electricity_generated - consumed_electricity
    return consumed_electricity
get_excess_electricity = lambda: excess_electricity[0]

metrics = (Metric('Internal rate of return', lc.lipidcane_tea.solve_IRR),
           Metric('Biodiesel production cost', get_biodiesel_prodcost, 'USD/yr'),
           Metric('Ethanol production cost', get_etoh_prodcost, 'USD/yr'),
           Metric('Fixed capital investment', get_FCI, 'USD'),
           Metric('Biodiesel production', get_biodiesel_prod, 'kg/hr'),
           Metric('Ethanol production', get_etoh_prod, 'kg/hr'),
           Metric('Steam', get_steam, 'MT/yr'),
           Metric('Consumed electricity', get_consumed_electricity, 'MWhr/yr'),
           Metric('Excess electricity', get_excess_electricity, 'MWhr/yr'))

lipidcane_model = Model(lc_sys, metrics)
lipidcane_model.load_default_parameters(lipid_cane)
param = lipidcane_model.parameter

# Lipid extraction rate
Mill = lc.system.pretreatment.U201
Mill_split = Mill.split
Lipid_index = Mill.outs[0].index('Lipid')
@param(element=Mill,
       distribution=triang(Mill_split[Lipid_index]),
       kind='coupled')
def set_lipid_extraction_rate(lipid_extraction_rate):
    Mill_split[Lipid_index] = lipid_extraction_rate
    
# Transesterification efficiency (both tanks)
R401 = lc.system.biodiesel.R401
@param(element=R401, distribution=triang(R401.efficiency), kind='coupled')
def set_transesterification_401_efficiency(efficiency):
    R401.efficiency = efficiency

R402 = lc.system.biodiesel.R402
@param(element=R402, distribution=triang(R402.efficiency), kind='coupled')
def set_transesterification_402_efficiency(efficiency):
    R402.efficiency = efficiency

# Fermentation efficiency
fermentation = lc.system.ethanol.R301
@param(element=fermentation, distribution=triang(fermentation.efficiency),
       kind='coupled')
def set_fermentation_efficiency(efficiency):
    fermentation.efficiency= efficiency
    
# Boiler efficiency
BT = lc.system.biorefinery.BT
@param(element=BT, distribution=triang(BT.boiler_efficiency))
def set_boiler_efficiency(boiler_efficiency):
    BT.boiler_efficiency = boiler_efficiency

# Turbogenerator efficiency
@param(element=BT, distribution=triang(BT.turbogenerator_efficiency))
def set_turbogenerator_efficiency(turbo_generator_efficiency):
    BT.turbo_generator_efficiency = turbo_generator_efficiency
    
# RVF separation
rvf = lc.system.pretreatment.C202
@param(element=rvf, distribution=triang(rvf.split['Lignin']),
        kind='coupled')
def set_rvf_solids_retention(solids_retention):
    rvf.split['Lignin', 'CaO', 'Ash', 'Cellulose', 'Hemicellulose'] = solids_retention

lipidcane_model_with_lipidfraction_parameter = lipidcane_model.copy()
lipidcane_model_with_lipidfraction_parameter.parameter(lc.set_lipid_fraction,
                                                       element=lipid_cane,
                                                       name='Lipid fraction',
                                                       distribution=triang(0.05))











