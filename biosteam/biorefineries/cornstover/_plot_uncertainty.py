# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import matplotlib.pyplot as plt
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_single_points, \
                                                 plot_montecarlo

    
# %% Plot Monte Carlo

light_color = colors.blue_tint.RGBn
dark_color = colors.blue_shade.RGBn
dot_color = colors.purple_shade.RGBn
position = 0
data = pd.read_excel('Monte Carlo cornstover.xlsx')
MESP = data['Minimum ethanol selling price'][2:]
bx = plot_montecarlo(MESP, light_color, dark_color, position)
sc = plot_single_points((position,), (2.15,), dot_color)

plt.ylim(1.65, 2.65)
plt.xlim(-1, 1)
plt.xticks([0], ['Cornstover'])
plt.ylabel('MESP ' + '[$\mathrm{USD} \cdot \mathrm{gal}^{-1}$]')
plt.tight_layout()
plt.legend([bx, sc], ['BioSTEAM', 'ASPEN (Humbird 2011)'])
