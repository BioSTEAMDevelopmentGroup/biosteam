# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_single_points, plot_horizontal_line, \
                                                 plot_montecarlo, plot_vertical_line

# light_color = colors.blue_tint.RGBn
# dark_color = colors.blue_shade.RGBn
# dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))

fig, axes = plt.subplots(ncols=1, nrows=2, constrained_layout=True, gridspec_kw=dict(height_ratios=[4, 1]))
metric_over_baseline_ax, relative_baseline_magnitude_ax = axes
plt.sca(metric_over_baseline_ax)

# %% Plot MESP and ethanol sales

positions_economic = (16, 17)
data = pd.read_excel('Monte Carlo cornstover.xlsx', header=[0, 1])
economic_data = np.array(data[[('Biorefinery', 'Minimum ethanol selling price'),
                               ('Biorefinery', 'Ethanol sales')]])
economic_data[:, 0] *= 100/2.15
economic_data[:, 1] *= 100/131e6
bx_economic = plot_montecarlo(economic_data,
                              colors.blue_tint.RGBn,
                              colors.blue_shade.RGBn, positions_economic, transpose=False)

# %% Plot electricity

plot_vertical_line(15.5, color=colors.grey_tint.RGBn)
positions_electricity = tuple(range(9))
areas = [f'Area {i}00' for i in nums]
units = 'MW'
positions = np.arange(0, 9)
electricity_cols = [(i, 'Electricity') for i in areas]
humbird_electricity = 41 * np.array((0.02, 0.14, 0.06, 0.05,
                                     0.18, 0.003, 0.03, 0.08,
                                     0.44))
electricity_data = data[electricity_cols]  #/ humbird_electricity
# electricity_data[('Consumed', 'Electricity')] = electricity_data.sum(1)
electricity_data[('Biorefinery', 'Excess electricity')] = data[('Biorefinery', 'Excess electricity')]/1000
electricity_data_humbird_normalized = electricity_data * (100/humbird_electricity)
bx_electricity = plot_montecarlo(electricity_data_humbird_normalized,
                                 colors.orange_tint.RGBn,
                                 colors.orange_shade.RGBn,
                                 transpose=True, positions=positions_electricity)

# plot_vertical_line(7.5, color=colors.orange_tint.shade(15).RGBn, ls='-.')
# plot_vertical_line(9.5, color=colors.orange_tint.shade(15).RGBn, ls='-.')

# %% Plot installation cost

units = '10^6 USD'
plot_vertical_line(8.5, color=colors.grey_tint.RGBn)
positions_installation = tuple(range(9, 16))
installation_cols = [(i, 'Installation cost') for i in areas[1:]]
humbird_installation = np.array([24.2, 32.9, 31.2, 22.3, 49.4, 5, 66, 6.9])
installation_data = data[installation_cols]
# installation_data[('Biorefinery', 'Installation cost')] = installation_data.sum(1)
installation_data_humbird_normalized = installation_data * (100/humbird_installation[1:]/1e6)
bx_installation = plot_montecarlo(installation_data_humbird_normalized,
                                  colors.purple_tint.RGBn, colors.purple_shade.RGBn,
                                  transpose=True, positions=positions_installation)
plot_horizontal_line(100, ls='--')
u_lb = 0; y_ub = 250
plt.ylim(0, 250)

y_text = 0.90*y_ub
plt.text(4, y_text, "Electricity demand", color=colors.orange_shade.RGBn,
         horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(12, y_text, "Installation cost", color=colors.purple_shade.RGBn,
         horizontalalignment='center', fontsize=12, fontweight='bold')
# plt.text(16, y_text, "Ethanol\nsales", color=colors.red_shade.RGBn,
#           horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(16.5, 0.975*y_text, "Economic\nindicators", color=colors.blue_shade.RGBn,
          horizontalalignment='center', fontsize=12, fontweight='bold')

plt.xlim(-0.5, 17.5)

plt.ylabel("Metric over baseline (%)")

plt.sca(relative_baseline_magnitude_ax)
plt.fill_between([-0.5, 15.5], 0, 1, color=colors.neutral_tint.tint(80).RGBn)

electricity_areas = humbird_electricity.copy()
electricity_areas /= max(electricity_areas)
plt.bar(positions_electricity, electricity_areas, 0.5,
        align='center', label="Electricity demand",
        color=colors.orange.tint(30).shade(15).RGBn,
        edgecolor=colors.orange_shade.RGBn)
plot_vertical_line(8.5, color=colors.grey_tint.RGBn)

installation_areas = humbird_installation[1:].copy()
installation_areas /= max(installation_areas)
plt.bar(positions_installation, installation_areas, 0.5,
        align='center', label="Installation cost",
        color=colors.purple.tint(30).shade(15).RGBn,
        edgecolor=colors.purple_shade.RGBn)

plot_vertical_line(15.5, color=colors.grey_tint.RGBn)
# plot_vertical_line(15.5, color='k', lw=0.8)
# plt.bar(positions_MESP, [1], 0.5,
#         align='center', label="MESP",
#         color=colors.blue.tint(30).shade(15).RGBn,
#         edgecolor=colors.blue_shade.RGBn)

plot_vertical_line(-0.5, color='k')
plt.hlines([0], [-0.5], [15.5], color='k')
relative_baseline_magnitude_ax.spines['bottom'].set_visible(False)
relative_baseline_magnitude_ax.spines['right'].set_visible(False)
plt.yticks([], [])
plt.ylim(0, 1)
plt.xlim(-0.5, 17.5)
area_marks = [i.replace(' ', '\n') for i in areas]
xmarks = area_marks + ['Excess\nProduced'] + area_marks[1:] 
plt.xticks(positions_electricity + positions_installation, xmarks)
plt.text(-0.25, 0.70, "Relative baseline magnitude", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=12, fontweight='bold')
plt.subplots_adjust(hspace=.0)

plt.sca(metric_over_baseline_ax)
metric_over_baseline_ax.set_zorder(1e6)
plt.xticks(positions_economic, ['Ethanol\nsales', 'MESP'])

# plot_vertical_line(15.5, color=colors.purple.tint(20).shade(10).RGBn, ls='-.')

# leg1 = ax.legend([bx_economic['boxes'][0]], ['MESP'], loc="upper left")
# leg2 = ax.legend([bx_electricity['boxes'][0]], ['Electricity'], loc="upper center")
# leg3 = ax.legend([bx_installation['boxes'][0]], ['Installation'], loc="upper right")
# ax.add_artist(leg2)
# ax.add_artist(leg1)

# light_box = Patch(color=colors.neutral_tint.RGBn)
# line = Line2D([0], [0], color=colors.neutral_shade.RGBn)
# dark_box = Patch(color=colors.neutral_shade.RGBn)
# plt.legend([(light_box, line), dark_box], ["Metric over baseline (%)", "Relative baseline magnintude"])

# leg1 = ax.legend([bx_economic['boxes'][0]], ['MESP'], loc="upper left")