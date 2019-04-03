# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 18:46:30 2019

@author: yoelr
"""

from biosteam import results_table, save_system_results
from lipidcane import lipidcane_sys

lipidcane_sys.simulate()
x = save_system_results(lipidcane_sys)
#x = results_table(lipidcane_sys.units)