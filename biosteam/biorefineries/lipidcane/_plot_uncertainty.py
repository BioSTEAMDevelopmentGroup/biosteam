# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 01:28:21 2019

@author: yoelr
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam.utils import DoubleColorCircle
from biosteam.evaluation.evaluation_tools import plot_single_points, \
     plot_montecarlo_across_coordinate, plot_montecarlo, plot_vertical_line, annotate_line
from biosteam import colors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# %% Constants
     
rho_etoh = 0.789 # kg/l
rho_bd = 0.870 # kg/l

# %% Plot montecarlo of lipidcane across lipid fraction

x_vals = [*range(1,11)]
x_labels = ['Sugar\ncane'] + x_vals
x_ticks = [0] + x_vals
def set_x_axis(with_labels=True):
    plt.xlim(-1, 10.5)
    if with_labels:
        plt.xticks(x_ticks, x_labels)
    else:
        plt.xticks(x_ticks, ())
        
# Plot metrics across lipid fraction

readxl = lambda sheet: pd.read_excel('Monte Carlo across lipid fraction.xlsx',
                                     sheet_name=sheet, index_col=0)

fig = plt.figure()

# IRR
IRR_ax = plt.subplot(3, 2, 1)
IRR = readxl('Internal rate of return') * 100 # To percent
lipid_fraction = np.array(IRR.columns) * 100
plt.ylabel('Internal rate of return [%]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, IRR)[2] # p50
annotate_line('IRR', 3, lipid_fraction, ys,
              dy=6, dy_text=0.8, position='over')

# TCI
TCI_ax = plt.subplot(3, 2, 3)
TCI = readxl('Fixed capital investment') * 1.05 / 1e6 # Account for working capital
plt.ylabel('Total capital investment [$10^6 \cdot \mathrm{USD}$]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, TCI)[2]
annotate_line('TCI', 3, lipid_fraction, ys, 
              dy=25, dy_text=2, position='over')

# Production
production_ax = plt.subplot(3, 2, 2)
ethanol_production = readxl('Ethanol production') / (1e6*rho_etoh)
biodiesel_production = readxl('Biodiesel production') / (1e6*rho_bd)
plt.ylabel('Production [$10^6 \cdot \mathrm{litter} \cdot \mathrm{yr}^{-1}$]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, ethanol_production,
                                       colors.orange_tint.RGBn,
                                       colors.orange_shade.RGBn)[2]
annotate_line('Ethanol', 8, lipid_fraction, ys, 
              dy=35, dy_text=5, position='over', color=colors.orange_shade.RGBn)
ys = plot_montecarlo_across_coordinate(lipid_fraction, biodiesel_production,
                                       colors.blue_tint.RGBn,
                                       colors.blue_shade.RGBn)[2]
annotate_line('Biodiesel', 4, lipid_fraction, ys, 
              dy=35, dy_text=4, position='over', color=colors.blue_shade.RGBn)

# Production cost
production_cost_ax = plt.subplot(3, 2, 4)
ethanol_production_cost = readxl('Ethanol production cost') / ethanol_production / 1e6
biodiesel_production_cost = readxl('Biodiesel production cost') / biodiesel_production / 1e6
plt.ylabel('Production cost [$\mathrm{USD} \cdot \mathrm{liter}^{-1}$]')
plot_montecarlo_across_coordinate(lipid_fraction, ethanol_production_cost,
                                  colors.orange_tint.RGBn,
                                  colors.orange_shade.RGBn)
plot_montecarlo_across_coordinate(lipid_fraction, biodiesel_production_cost,
                                  colors.blue_tint.RGBn,
                                  colors.blue_shade.RGBn)

# Steam
steam_ax = plt.subplot(3, 2, 5)
steam = readxl('Steam')/1000
plt.ylabel('Steam [$10^{3} \cdot \mathrm{MT} \cdot \mathrm{yr}^{-1}$]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, steam)[2]
annotate_line('Steam', 8, lipid_fraction, ys, 
              dy=150, dy_text=20, position='over')


# Electricity
electricity_ax = plt.subplot(3, 2, 6)
consumed_electricity = readxl('Consumed electricity')/1000
excess_electricity = readxl('Excess electricity')/1000
plt.ylabel('Electricity [$\mathrm{GWhr} \cdot \mathrm{yr}^{-1}$]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, consumed_electricity,
                                       colors.purple_tint.RGBn,
                                       colors.purple_shade.RGBn)[2]
annotate_line('Consumed electricity', 7.5, lipid_fraction, ys, 
              dy=90, dy_text=20, position='over',
              color=colors.purple_shade.RGBn)
ys = plot_montecarlo_across_coordinate(lipid_fraction, excess_electricity,
                                       colors.yellow_tint.RGBn,
                                       colors.yellow_shade.RGBn)[2]
annotate_line('Excess electricity', 4, lipid_fraction, ys, 
              dy=150, dy_text=20, position='over',
              color=colors.yellow_shade.RGBn)


# Plot sugarcane values and SuperPro values
x_superpro = [0, 2, 5, 10]
data_sc = pd.read_excel('Monte Carlo sugarcane.xlsx', header=[0,1])
get_metric = lambda name: np.asarray(data_sc['Biorefinery', name]).flatten()

# IRR
plt.sca(IRR_ax)
plot_single_points(x_superpro, [13.5, 13.7, 15.2, 17.5])
IRR = get_metric('Internal rate of return') * 100 # To percent
plot_montecarlo(IRR)
plot_vertical_line(1)
IRR_ub = 30
plt.ylim(0, IRR_ub)
IRR_yticks = np.arange(0, 31, 30/5)
plt.yticks(IRR_yticks)
set_x_axis(False)
y_text = 0.85*IRR_ub
plt.text(0.05, y_text, "A", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')


# TCI
plt.sca(TCI_ax)
TCI = get_metric('Fixed capital investment')  * 1.05 / 1e6 # Account for working capital
plot_single_points(x_superpro, [158.5, 172.9, 178.3, 195.0])
plot_montecarlo(TCI)
plot_vertical_line(1)
TCI_ub = 300
plt.ylim(0, TCI_ub)
TCI_yticks = np.arange(0, 300, 300/5)
plt.yticks(TCI_yticks)
set_x_axis(False)
y_text = 0.85*TCI_ub
plt.text(0.05, y_text, "B", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')


# Production
plt.sca(production_ax)
plot_single_points(x_superpro[-1], [48], colors.blue_shade.RGBn)
plot_single_points([0, x_superpro[-1]], [141, 70], colors.orange_shade.RGBn)
ethanol_production = get_metric('Ethanol production') / (1e6*rho_etoh)
plot_montecarlo(ethanol_production,
                colors.orange_tint.RGBn,
                colors.orange_shade.RGBn)
production_ub = 225
plt.ylim(0, production_ub)
production_yticks = np.arange(0, 226, 225/5)
plt.yticks(production_yticks)
set_x_axis(False)
plot_vertical_line(1)
y_text = 0.85*production_ub
plt.text(0.05, y_text, "D", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')


# Production cost
plt.sca(production_cost_ax)
plot_single_points(x_superpro[1:], [0.89, 0.84, 0.76], colors.blue_shade.RGBn)
plot_single_points(x_superpro, [0.48, 0.46, 0.44, 0.4], colors.orange_shade.RGBn)
ethanol_production_cost = get_metric('Ethanol production cost') / ethanol_production / 1e6
plot_montecarlo(ethanol_production_cost,
                colors.orange_tint.RGBn,
                colors.orange_shade.RGBn)
plot_vertical_line(1)
production_cost_ub = 1.2
plt.ylim(0, production_cost_ub)
production_cost_yticks = np.arange(0, 1.2, 1.2/5)
plt.yticks(production_cost_yticks)
set_x_axis(False)
y_text = 0.85*production_cost_ub
plt.text(0.05, y_text, "E", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')


# Steam
plt.sca(steam_ax)
plot_single_points([0, 10], [686.056, 656.000])
steam = get_metric('Steam')/1000
plot_montecarlo(steam)
plot_vertical_line(1)
steam_ub = 900
plt.ylim(0, steam_ub)
steam_yticks = np.arange(0, 900, 900/5)
plt.yticks(steam_yticks)
set_x_axis(True)
plt.xlabel('Feedstock lipid content [%]')
y_text = 0.85*steam_ub
plt.text(0.05, y_text, "C", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')


# Electricity
plt.sca(electricity_ax)
plot_single_points([0, 10], [50.187, 62.644], colors.purple_shade.RGBn)
plot_single_points([0, 10], [110, 260], colors.yellow_shade.RGBn)
consumed_electricity = get_metric('Consumed electricity')/1000
excess_electricity = get_metric('Excess electricity')/1000
plot_montecarlo(consumed_electricity,
                colors.purple_tint.RGBn,
                colors.purple_shade.RGBn)
plot_montecarlo(excess_electricity,
                colors.yellow_tint.RGBn,
                colors.yellow_shade.RGBn)
plot_vertical_line(1)
electricity_ub = 575
plt.ylim(0, electricity_ub)
electricity_yticks = np.arange(0, 575, 575/5)
plt.yticks(electricity_yticks)
set_x_axis(True)
plt.xlabel('Feedstock lipid content [%]')
y_text = 0.85*electricity_ub
plt.text(0.05, y_text, "F", color=colors.neutral_shade.RGBn,
         horizontalalignment='center', fontsize=14, fontweight='bold')



plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=0.3)
IRR_ax.tick_params(axis="x", direction="inout", length=4)
IRR_ax.set_zorder(2)
TCI_ax.tick_params(axis="x", direction="inout", length=4)
TCI_ax.set_zorder(1)
production_ax.tick_params(axis="x", direction="inout", length=4)
production_ax.set_zorder(2)
production_cost_ax.tick_params(axis="x", direction="inout", length=4)
production_cost_ax.set_zorder(1)
steam_ax.tick_params(axis="x", direction="inout", length=4)
electricity_ax.tick_params(axis="x", direction="inout", length=4)
axs = [TCI_ax, IRR_ax, production_ax, production_cost_ax, electricity_ax, steam_ax]
yticks = [TCI_yticks, IRR_yticks, production_yticks, 
          production_cost_yticks, electricity_yticks, steam_yticks]
ubs = [TCI_ub, IRR_ub, production_ub, 
       production_cost_ub, electricity_ub, steam_ub]
ax2s = [ax.twinx() for ax in axs]
for ub, yt, ax, ax2 in zip(ubs, yticks, axs, ax2s):
    ax.tick_params(axis='y', right=False, direction="inout", length=4)
    plt.sca(ax2)
    plt.yticks(yt, ['']*len(yt))
    plt.ylim(0, ub)
    ax2.zorder = 1000
    ax2.tick_params(direction="in")

for ax in (IRR_ax, production_ax):
    ax2 = ax.twiny()
    plt.sca(ax2)
    plt.xticks(x_ticks, ())
    ax2.zorder = 1000
    ax2.tick_params(direction="in")


plt.sca(production_ax)
BioSTEAM_patch = Patch(facecolor=colors.neutral_tint.RGBn, 
                       edgecolor=colors.neutral_shade.RGBn,
                       label='BioSTEAM')
Baseline_circle = Line2D([0], [0], marker='o',
                         color='w',
                         label='Baseline',
                         markerfacecolor=colors.neutral_shade.RGBn,
                         markersize=10)
legend = production_ax.legend(handles=[BioSTEAM_patch, Baseline_circle])
frame = legend.get_frame()
frame.set_linewidth(0.0)
frame.set_facecolor('none')

fig.align_ylabels([IRR_ax, TCI_ax, steam_ax])
fig.align_ylabels([production_ax, production_cost_ax, electricity_ax])
fig.align_ylabels([IRR_ax, TCI_ax])

# plt.sca(electricity_ax)
# legend = DoubleColorLegend()
# legend.add_box('BioSTEAM',
#                leftcolor=colors.purple_tint.RGBn, 
#                rightcolor=colors.yellow_tint.RGBn)
# legend.add_circle('SuperPro (Huang 2016)',
#                   leftcolor=colors.purple_shade.RGBn, 
#                   rightcolor=colors.yellow_shade.RGBn)
# legend.legend()