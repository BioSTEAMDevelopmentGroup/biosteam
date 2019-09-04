# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:15 2019

@author: yoelr
"""
import biosteam as bst

bst.CE = 525

# kg/hr
factor = 1/907.18474
ethanol_density_kgL = 0.789 # kg/L
liter_per_gallon = 3.78541
ethanol_cost = 2.15 # USD/gal
ethanol_density_kggal = liter_per_gallon * ethanol_density_kgL # kg/gal

price = {'Ethanol': ethanol_cost/ethanol_density_kggal,
         'Feedstock': 46.8 * factor,
         'Sulfuric acid': 81.39 * factor,
         'Ammonia': 406.96 * factor,
         'CSL': 51.55 * factor,
         'DAP': 895.32 * factor,
         'Sorbitol': 1021.93 * factor,
         'Glucose': 526.52 * factor,
         'Caustic': 135.65 * factor,
         'Boiler chems': 4532.17 * factor,
         'FGD lime': 180.87 * factor,
         'Cooling tower chems': 2716.1 * factor,
         'Makeup water': 0.23 * factor,
         'Ash disposal': -28.86 * factor,
         'Electricity': 0.0572, # USD/kWh
         'Denaturant': 0.756,
         'Cellulase': 0.21} 
bst.PowerUtility.price = price['Electricity']
_ha = bst.HeatUtility.heating_agents['Low pressure steam']
_ha.efficiency = 0.80
_ha.T = 529.2
_ha.P = 44e5
_ha.price_kmol = 0.30626
_ha.Hvap = 30235.86
_CW = bst.HeatUtility.cooling_agents['Cooling water']
_CW.T = 28 + 273.15
_CW.T_limit = _CW.T + 9
_CW.price_kmol = 0
bst.HeatUtility.cooling_agents['Chilled water'].price_kJ = 0
