# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:11:08 2019

@author: yoelr
"""

from biosteam import PowerUtility, HeatUtility, find, Flowsheet
import biosteam

__all__ = ('price',)

# %% Process settings

biosteam.CE = 567 # 2013
find.set_flowsheet(Flowsheet('sugarcane'))
PowerUtility.price = 0.065
_ha = HeatUtility.heating_agents['Low pressure steam']
_ha.efficiency = 0.90
_ha.T = 529.2
_ha.P = 44e5
_ha.price_kmol = 0.30626
_ha.Hvap = 30235.86
HeatUtility.cooling_agents['Cooling water'].price_kmol = 0
HeatUtility.cooling_agents['Chilled water'].price_kJ = 0


# Raw material price (USD/kg)
price = {'Sugar cane': 0.03455, # 70% m.c
         'Water': 0.000353,
         'HCl': 0.205,
         'Lime': 0.077,
         'H3PO4': 0,#0.700, # TODO: find price
         'NaOH':0.41,
         'Protease': 0.5,
         'Polymer': 0, # TODO: find price
         'Steam': 0.017,
         'Ethanol': 0.789,
         'Waste': -0.33,
         'Gasoline': 0.756} # 2 USD/gal