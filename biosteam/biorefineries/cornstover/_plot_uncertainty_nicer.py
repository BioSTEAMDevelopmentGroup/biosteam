# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_single_points, plot_horizontal_line, \
                                                 plot_montecarlo, plot_vertial_line

light_color = colors.blue_tint.RGBn
dark_color = colors.blue_shade.RGBn
dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))
    
# %% Plot MESP

plt.figure()

positions_economic = (0,)
data = pd.read_excel('Monte Carlo cornstover.xlsx', header=[0, 1])
MESP = np.array(data[('Biorefinery', 'Minimum ethanol selling price')])*(100/2.15)
bx_economic = plot_montecarlo(MESP, light_color, dark_color, positions_economic, transpose=True)

# %% Plot electricity

plot_vertial_line(0.5)
positions_electricity = tuple(range(1, 10))
areas = [f'Area {i}00' for i in nums]
units = 'MW'
positions = np.arange(0, 9)
electricity_cols = ([(i, 'Electricity') for i in areas]
                    + [('Biorefinery', 'Excess electricity')])
humbird_electricity = 41 * np.array((0.02, 0.14, 0.06, 0.05,
                                     0.18, 0.003, 0.03, 0.08, 0.44))
electricity_data = data[electricity_cols]  #/ humbird_electricity
electricity_data[('Biorefinery', 'Excess electricity')] /= 1000
electricity_data = electricity_data * (100/humbird_electricity)
bx_electricity = plot_montecarlo(electricity_data,
                                 colors.orange_tint.RGBn,
                                 colors.orange_shade.RGBn,
                                 transpose=True, positions=positions_electricity)
area_marks = [i.replace(' ', '\n') for i in areas]
xmarks = ['Biorefinery'] + area_marks + ['Excess'] + area_marks


# %% Plot installation cost

units = '10^6 USD'
plot_vertial_line(9.5)
positions_installation = tuple(range(10, 17))
installation_cols = [(i, 'Installation cost') for i in areas[1:]]
humbird_installation = np.array([24.2, 32.9, 31.2, 22.3, 49.4, 5, 66, 6.9])
installation_data = data[installation_cols] * (100/humbird_installation[1:]/1e6)
bx_installation = plot_montecarlo(installation_data, colors.purple_tint.RGBn,
                                  colors.purple_shade.RGBn, transpose=True,
                                  positions=positions_installation)
# plot_single_points(positions[:-2], humbird_installation[1:], dot_color)
plot_horizontal_line(100, ls='--')
plt.ylim(0, 250)
plt.xticks(positions_economic + positions_electricity + positions_installation, xmarks)
def m(j):
    try:
        return j[0]
    except:
        return j
bx_electricity = {i:m(j) for i,j in bx_electricity.items()}
bx_installation = {i:m(j) for i,j in bx_installation.items()}
plt.legend([bx_economic, bx_electricity, bx_installation], ['MESP', 'Electricity', 'Installation'])
plt.xlim(-0.5, 17)