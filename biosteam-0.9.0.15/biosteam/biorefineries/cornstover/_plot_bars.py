# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_single_points, plot_bars
from biosteam.biorefineries.cornstover.model import metrics

light_color = colors.blue_tint.RGBn
dark_color = colors.blue_shade.RGBn
dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))
data = pd.Series([i() for i in metrics], index=[i.index for i in metrics])


# %% Plot electricity

areas = [f'Area {i}00' for i in nums]
xmarks = [i.replace(' ', '\n') for i in areas]
units = 'MW'
electricity_cols = [(i, 'Electricity') for i in areas]
humbird_electricity = 41 * np.array((0.02, 0.14, 0.06, 0.05,
                                     0.18, 0.003, 0.03, 0.08))
electricity_data = data[electricity_cols]  #/ humbird_electricity

labels = ('BioSTEAM', 'Humbird (Aspen)')
ys = (electricity_data, humbird_electricity)
colors = (dark_color, dot_color)
edgecolors = ('k', 'k')
scenarios = areas
plt.figure()
plot_bars(scenarios, ys, colors,
          edgecolors, labels)
plt.ylabel('Electricity (MW)')

# %% Plot installation cost

units = '10^6 USD'
installation_cols = [(i, 'Installation cost') for i in areas]
humbird_installation = np.array([24.2, 32.9, 31.2, 22.3, 49.4, 5, 66, 6.9])
installation_data = data[installation_cols] / 1e6

labels = ('BioSTEAM', 'Humbird (Aspen)')
ys = (installation_data[1:],
      humbird_installation[1:])
plt.figure()
plot_bars(scenarios[1:], ys, colors,
          edgecolors, labels)
plt.ylabel('Installation cost ($10^6$ USD)')



# # %% Plot cooling duty

# units = 'MMkcal/hr'
# cooling_cols = [(i, 'Cooling duty') for i in areas]
# cooling_data = data['Cooling duty']
# humbird_cooling = 87 * np.array([])
# ys = (humbird_cooling,
#       cooling_data)
# plt.figure()
# plot_bars(scenarios, ys, colors,
#           edgecolors, labels)


