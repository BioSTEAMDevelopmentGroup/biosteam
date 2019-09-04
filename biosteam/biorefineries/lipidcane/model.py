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
Ethanol = lc.system.ethanol.Ethanol
Biodiesel = lc.system.biodiesel.Biodiesel
Lipid_cane = lc.system.pretreatment.Lipid_cane
etoh_prodcost = [0]
def get_biodiesel_prodcost():
    bd, etoh_prodcost[0] = tea.production_cost(Biodiesel, Ethanol)
    return bd
get_etoh_prodcost = lambda: etoh_prodcost[0]
get_FCI = lambda: tea._FCI_cached
etoh_prod = [0]
def get_biodiesel_prod():
    bd, etoh_prod[0] = np.array([Biodiesel.massnet, Ethanol.massnet]) * tea._annual_factor
    return bd
get_etoh_prod = lambda: etoh_prod[0]

metrics = (Metric('Internal rate of return', '', lc.lipidcane_tea.solve_IRR),
           Metric('Biodiesel production cost', 'USD/yr', get_biodiesel_prodcost),
           Metric('Ethanol production cost', 'USD/yr', get_etoh_prodcost),
           Metric('Fixed capital investment', 'USD', get_FCI),
           Metric('Biodiesel production', 'kg/hr', get_biodiesel_prod),
           Metric('Ethanol production', 'kg/hr', get_etoh_prod))

lipidcane_model = Model(lc.lipidcane_sys, metrics)
lipidcane_model.load_default_parameters(Lipid_cane)
param = lipidcane_model.parameter

# Lipid extraction rate
Mill = lc.system.pretreatment.Mill
Mill_split = Mill.split
Lipid_index = Mill.outs[0].index('Lipid')
@param(element=Mill,
       distribution=triang(Mill_split[Lipid_index]),
       kind='coupled')
def set_lipid_extraction_rate(lipid_extraction_rate):
    Mill_split[Lipid_index] = lipid_extraction_rate
    
# Transesterification efficiency (both tanks)
trans_reactors = [lc.system.biodiesel.R1, lc.system.biodiesel.R2]
for Ri in trans_reactors:
    @param(element=Ri, distribution=triang(Ri.efficiency), kind='coupled')
    def set_transesterification_efficiency(efficiency):
        for i in trans_reactors:
            i.efficiency = efficiency

# Fermentation efficiency
fermentation = lc.system.ethanol.P24
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
rvf = lc.system.pretreatment.P14
@param(element=rvf, distribution=triang(rvf.split['Lignin']),
        kind='coupled')
def set_rvf_solids_retention(solids_retention):
    rvf.split['Lignin', 'CaO', 'Ash', 'Cellulose', 'Hemicellulose'] = solids_retention

lipidcane_model_with_lipidfraction_parameter = lipidcane_model.copy()
lipidcane_model_with_lipidfraction_parameter.parameter(lc.set_lipid_fraction,
                                                       element=Lipid_cane,
                                                       name='Lipid fraction',
                                                       distribution=triang(0.05))











