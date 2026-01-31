# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 10:54:34 2025

@author: yoelr
"""
from pint import UnitRegistry
u = UnitRegistry()
yogurt_wt = 4.9 * u.lb
yogurt_density = 1.03 * u.g / u.mL
print((yogurt_wt / yogurt_density).to('gal'))
