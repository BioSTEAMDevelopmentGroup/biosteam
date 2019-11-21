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

data = pd.read_excel('Monte Carlo cornstover.xlsx', header=[0, 1])

# %% Plot MESP
# plt.figure()

# posistions_MESP = (0,)
# MESP_data = data[('Biorefinery', 'Minimum ethanol selling price')]
# bx_MESP = plot_montecarlo(MESP_data,
#                           colors.blue_tint.RGBn,
#                           colors.blue_shade.RGBn, 
#                           posistions_MESP, transpose=False)
# dot_MESP = plot_single_points(posistions_MESP, [2.15], s=125, color=colors.red_dark.RGBn)
# plt.ylabel('MESP ($\mathrm{USD} \cdot \mathrm{gal}^{-1}$)')
# plt.ylim(1.5, 2.50)
# bx_patch = Patch(facecolor=colors.blue_tint.RGBn, edgecolor=colors.blue_shade.RGBn)
# plt.legend([bx_patch, dot_MESP], ['BioSTEAM', 'benchmark'])
# plt.xticks([], [])
                                                 
# %% Setup of subplots

# light_color = colors.blue_tint.RGBn
# dark_color = colors.blue_shade.RGBn
# dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))

fig, axes = plt.subplots(ncols=1, nrows=2, constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 4]))
magnitude_ax, metric_over_benchmark_ax = axes
plt.sca(metric_over_benchmark_ax)

# %% Plot MESP and ethanol sales

positions_other = (16, 17, 18)
other_index = [('Biorefinery', 'Steam demand'),
               ('Biorefinery', 'Ethanol production'),
               ('Biorefinery', 'Minimum ethanol selling price')]
other_data = np.array(data[other_index])
other_data[:, 0] *= 100./234784.
other_data[:, 1] *= 100./22273.
other_data[:, 2] *= 100./2.15
bx_other = plot_montecarlo(other_data,
                           colors.blue_tint.RGBn,
                           colors.blue_shade.RGBn, positions_other, transpose=False)

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
yticks = np.arange(0, 251, 50)
plt.yticks(yticks)
y_text = 0.885*y_ub
y_letter = 0.875 * y_ub
plt.text(4, y_text, "Electricity demand", color=colors.orange_shade.RGBn,
          horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(-0.25, y_letter, "C", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(12, y_text, "Installation cost", color=colors.purple_shade.RGBn,
          horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(8.75, y_letter, "D", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(15.75, y_letter, "E", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
# plt.text(9.25, y_letter, "D", color=colors.neutral_shade.RGBn,
#          horizontalalignment='center', fontsize=14, fontweight='bold')

# plt.text(16, y_text, "Ethanol\nsales", color=colors.red_shade.RGBn,
#           horizontalalignment='center', fontsize=12, fontweight='bold')
# plt.text(16.0, y_text, "MESP", color=colors.blue_shade.RGBn,
#           horizontalalignment='center', fontsize=12, fontweight='bold')

plt.xlim(-0.5, 18.5)
plt.ylabel("Metric over benchmark [%]")
area_marks = [i.replace(' ', '\n') for i in areas]
xmarks = area_marks + ['Excess'] + area_marks[1:] + ['Steam\ndemand', '   EtOH\n    prod.', '  MESP']
xticks = positions_electricity + positions_installation + positions_other
plt.xticks(xticks, xmarks)
metric_over_benchmark_ax.set_zorder(1e6)

plt.sca(magnitude_ax)
plt.fill_between([-0.5, 15.5], 0, 1, color=colors.neutral_tint.tint(85).RGBn)

electricity_areas = humbird_electricity # electricity_data.median()
electricity_areas /= max(electricity_areas)
plt.bar(positions_electricity, electricity_areas, 0.5,
        align='center', label="Electricity demand",
        color=colors.orange.tint(30).shade(15).RGBn,
        edgecolor=colors.orange_shade.RGBn)
plot_vertical_line(8.5, color=colors.grey_tint.RGBn)

installation_areas = humbird_installation[1:]# installation_data.median()
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
plt.hlines([1], [-0.5], [15.5], color='k')
magnitude_ax.spines['top'].set_visible(False)
magnitude_ax.spines['right'].set_visible(False)
magnitude_ax.tick_params(axis="x", direction="in", length=0)
magnitude_ax.set_zorder(2)
metric_over_benchmark_ax.set_zorder(1)
plt.yticks([], [])
plt.xticks(xticks[:-3], [])
plt.ylim(0, 1)
plt.xlim(-0.5, 18.5)
plt.text(4, 0.65, "Benchmark magnitude", color=colors.neutral.shade(35).RGBn,
         horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(-0.25, 0.60, "A", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(8.75, 0.60, "B", color=colors.neutral_shade.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.subplots_adjust(hspace=.0)

metric_over_benchmark_ax.tick_params(axis='x', direction="inout", length=4)
for ax in axes:
    ax.tick_params(axis='y', right=False, direction="inout", length=4)
ax2 = metric_over_benchmark_ax.twinx()
plt.sca(ax2)
plt.ylim(0, 250)
plt.yticks(yticks, [])
ax2.zorder = 1000
ax2.tick_params(direction="in")
    
xlabels = metric_over_benchmark_ax.get_xticklabels()    
# for xtick in xlabels[-3:]:
#     xtick.set_rotation(90)


# plot_vertical_line(15.5, color=colors.purple.tint(20).shade(10).RGBn, ls='-.')

# leg1 = ax.legend([bx_economic['boxes'][0]], ['MESP'], loc="upper left")
# leg2 = ax.legend([bx_electricity['boxes'][0]], ['Electricity'], loc="upper center")
# leg3 = ax.legend([bx_installation['boxes'][0]], ['Installation'], loc="upper right")
# ax.add_artist(leg2)
# ax.add_artist(leg1)

# light_box = Patch(color=colors.neutral_tint.RGBn)
# line = Line2D([0], [0], color=colors.neutral_shade.RGBn)
# dark_box = Patch(color=colors.neutral_shade.RGBn)
# plt.legend([(light_box, line), dark_box], ["Metric over benchmark (%)", "Relative benchmark magnintude"])

# leg1 = ax.legend([bx_economic['boxes'][0]], ['MESP'], loc="upper left")