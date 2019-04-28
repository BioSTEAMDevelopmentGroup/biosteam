# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 12:43:00 2019

@author: yoelr
"""

# %% Parameter values

import biosteam as bs
import lipidcane as lc

from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np

problem = {
    'num_vars': 2,
    'names': ['annual_capacity', 'lipid_content'],
    'bounds': [[0.5, 1.5],
               [0.02, 0.15]]
}
annual_capacity, lipid_content = saltelli.sample(problem, 3).transpose()

# %% Set Parameters

ss = bs.Sensitivity(lc.lipidcane_sys, lc.lipidcane_tea.solve_IRR, ID='IRR')

Lipid_cane = lc.Lipid_cane
original_flow = np.array(Lipid_cane.mol)
def change_flow(annual_capacity):
    Lipid_cane.mol = original_flow*annual_capacity

ss.addparam(Lipid_cane,
            lc.set_lipid_fraction,
            lipid_content)

ss.addparam(Lipid_cane,
            change_flow,
            annual_capacity)

ss.simulate(False)
Si = sobol.analyze(problem, ss.table['IRR'])
print(Si)
