# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 10:02:05 2019

@author: yoelr
"""

from biosteam import PowerUtility, HeatUtility, find, Flowsheet
import biosteam

__all__ = ('price',)

# %% Process settings

biosteam.CE = 567 # 2013
find.set_flowsheet(Flowsheet('lipidcane'))
_ha = HeatUtility.heating_agents['Low pressure steam']
_ha.efficiency = 0.85
_ha.T = 529.2
_ha.P = 44e5
_ha.price_kmol = 0.30626
_ha.Hvap = 30235.86
HeatUtility.cooling_agents['Cooling water'].price_kmol = 0
HeatUtility.cooling_agents['Chilled water'].price_kJ = 0

# Raw material price (USD/kg)
price = {'Electricity': 0.065,
         'Lipid cane': 0.03455, # 70% m.c
         'Methanol': 0.547,
         'Water': 0.000353,
         'HCl': 0.205,
         'Lime': 0.077,
         'H3PO4': 0, # 0.700, # TODO: find price
         'NaOCH3': 2.93,
         'NaOH':0.41,
         'Protease': 0.5,
         'Polymer': 0, # TODO: find price
         'Steam': 0.017,
         'Natural gas': 0.218,
         'Crude glycerol': 0.21,
         'Biodiesel': 1.38, # assuming density = 870 kg/m3
         'Ethanol': 0.789,
         'Waste': -0.33,
         'Gasoline': 0.756} # 2 USD/gal

PowerUtility.price = price['Electricity']
                
                

