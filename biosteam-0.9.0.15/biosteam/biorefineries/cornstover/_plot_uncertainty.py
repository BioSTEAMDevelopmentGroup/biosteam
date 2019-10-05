# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_single_points, \
                                                 plot_montecarlo

light_color = colors.blue_tint.RGBn
dark_color = colors.blue_shade.RGBn
dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))
    
# %% Plot MESP

plt.figure()

positions = (0,)
data = pd.read_excel('Monte Carlo cornstover.xlsx', header=[0, 1])
MESP = np.array(data[('Biorefinery', 'Minimum ethanol selling price')])*(100/2.15)
bx = plot_montecarlo(MESP, light_color, dark_color, positions, transpose=True)

plt.xticks([0], ['Cornstover'])
plt.ylim(50, 150)
plt.ylabel('MESP ' + '[$\mathrm{USD} \cdot \mathrm{gal}^{-1}$]')
plt.tight_layout()
plt.legend([bx], ['BioSTEAM'])

# %% Plot electricity


areas = [f'Area {i}00' for i in nums]
xmarks = [i.replace(' ', '\n') for i in areas]
units = 'MW'
positions = np.arange(0, 9)
electricity_cols = ([(i, 'Electricity') for i in areas]
                    + [('Biorefinery', 'Excess electricity')])
humbird_electricity = 41 * np.array((0.02, 0.14, 0.06, 0.05,
                                     0.18, 0.003, 0.03, 0.08, 0.44))
electricity_data = data[electricity_cols]  #/ humbird_electricity
electricity_data[('Biorefinery', 'Excess electricity')] /= 1000
electricity_data = electricity_data * (100/humbird_electricity)
plot_montecarlo(electricity_data, light_color, dark_color, transpose=True)
plt.ylim(0, 250)
plt.xticks(positions, xmarks + ['Excess'])
# plt.xlim(-0.5, 7.5)

# %% Plot installation cost

units = '10^6 USD'
installation_cols = [(i, 'Installation cost') for i in areas[1:]]
humbird_installation = np.array([24.2, 32.9, 31.2, 22.3, 49.4, 5, 66, 6.9])
installation_data = data[installation_cols] * (100/humbird_installation[1:]/1e6)
plot_montecarlo(installation_data, light_color, dark_color, transpose=True)
# plot_single_points(positions[:-2], humbird_installation[1:], dot_color)
plt.xticks(positions[:-2], xmarks[1:])
plt.ylim(0, 250)
plt.xlim(-0.5, 6.5)