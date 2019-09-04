# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 01:28:21 2019

@author: yoelr
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam.utils import DoubleColorLegend
from biosteam.evaluation.evaluation_tools import plot_spearman, plot_single_points, \
     plot_montecarlo_across_coordinate, plot_montecarlo, plot_vertial_line, annotate_line
from biosteam import colors

# %% Constants
     
rho_etoh = 0.789 # kg/l
rho_bd = 0.870 # kg/l


# %%  Plot spearman correlations

# Replacement parameter labels
replacement_labels = {
        'Stream-Ethanol price': 'Ethanol price',
        'TEA operating days': 'Operating days',
        'Stream-Lipid cane price': 'Lipid cane price',
        'Fermentation-U34 efficiency': 'Fermentation efficiency',
        'Stream-Biodiesel price': 'Biodiesel price',
        'Stream-Crude glycerol price': 'Crude glycerol price',
        'Stream-Lime price': 'Lime price',
        'Stream-Lipid cane flow rate': 'Lipid cane flow rate',
        'TEA income tax': 'Income tax',
        'Boiler turbogenerator-BT boiler efficiency': 'Boiler efficiency',
        'Boiler turbogenerator boiler exponent': 'Boiler cost exponent',
        'Boiler turbogenerator boiler base cost': 'Boiler base cost',
        'Boiler turbogenerator turbogenerator electricity rate': 'Turbogenerator electricity rate',
        'Boiler turbogenerator turbogenerator base cost': 'Turbogenerator base cost',
        'Mix tank exponent': 'Mix tank cost exponent',
        'Stream-Lipid cane lipid fraction': 'Lipid cane lipid fraction',
        'Transesterification-U61 efficiency': 'Transesterification efficiency',
        'Power utility price': 'Electricity price',
        'Cooling tower base cost': 'Cooling tower base cost'}

def replace_label_text(label_text):
    """Replace label text for graph."""
    name, distribution = label_text.split(' (')
    lb, mid, ub = eval('(' + distribution)
    if 'efficiency' in name:
        distribution = f" ({lb:.2f}, {mid:.2f}, {ub:.2f})"
    else:
        distribution = f" ({lb:.3g}, {mid:.3g}, {ub:.3g})"
    pos = name.find(' [')
    if pos != -1:
        units = name[pos:]
        name = name[:pos]
    else:
        units = ''
    if name in replacement_labels:
        name = replacement_labels[name]
    return name + units + distribution

# Get data
rhos = pd.read_excel('Spearman correlation lipidcane.xlsx',
                     header=[0])['Internal rate of return']

# Sort most important
rhos = rhos[(rhos.abs()>0.1)] 

# Plot and fix axis labels
fig, ax = plot_spearman(rhos)
labels = [item.get_text() for item in ax.get_yticklabels()]
new_labels = [replace_label_text(i) for i in labels]
ax.set_yticklabels(new_labels)

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
                                     sheet_name=sheet)

IRR_ax = plt.subplot(2, 2, 1)
IRR = readxl('Internal rate of return') * 100 # To percent
lipid_fraction = np.array(IRR.columns) * 100
plt.ylabel('Internal Rate of Return [%]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, IRR)[2] # p50
annotate_line('IRR', 3, lipid_fraction, ys,
              dy=6, dy_text=0.8, position='over')

TCI_ax = plt.subplot(2, 2, 3)
TCI = readxl('Fixed capital investment') * 1.05 / 1e6 # Account for working capital
plt.ylabel('Total Capital Investment [$10^6 \cdot \mathrm{USD}$]')
ys = plot_montecarlo_across_coordinate(lipid_fraction, TCI)[2]
annotate_line('TCI', 3, lipid_fraction, ys, 
              dy=25, dy_text=2, position='over')

production_ax = plt.subplot(2, 2, 2)
ethanol_production = readxl('Ethanol production') / (1e6*rho_etoh)
biodiesel_production = readxl('Biodiesel production') / (1e6*rho_bd)
plt.ylabel('Production [$10^6 \cdot \mathrm{litter}$]')
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

production_cost_ax = plt.subplot(2, 2, 4)
ethanol_production_cost = readxl('Ethanol production cost') / ethanol_production / 1e6
biodiesel_production_cost = readxl('Biodiesel production cost') / biodiesel_production / 1e6
plt.ylabel('Production cost [$\mathrm{USD} \cdot \mathrm{liter}^{-1}$]')
plot_montecarlo_across_coordinate(lipid_fraction, ethanol_production_cost,
                                  colors.orange_tint.RGBn,
                                  colors.orange_shade.RGBn)
plot_montecarlo_across_coordinate(lipid_fraction, biodiesel_production_cost,
                                  colors.blue_tint.RGBn,
                                  colors.blue_shade.RGBn)

# Plot sugarcane values and SuperPro values
x_superpro = [0, 2, 5, 10]
data_sc = pd.read_excel('Monte Carlo sugarcane.xlsx', header=[0,1])
get_metric = lambda name: np.asarray(data_sc[name]).flatten()

plt.sca(IRR_ax)
plot_single_points(x_superpro, [13.5, 13.7, 15.2, 17.5])
IRR = get_metric('Internal rate of return') * 100 # To percent
plot_montecarlo(IRR, position=0)
plot_vertial_line(1)
plt.ylim(0, 30)
set_x_axis(False)

plt.sca(TCI_ax)
TCI = get_metric('Fixed capital investment')  * 1.05 / 1e6 # Account for working capital
plot_single_points(x_superpro, [158.5, 172.9, 178.3, 195.0])
plot_montecarlo(TCI, position=0)
plot_vertial_line(1)
plt.ylim(120, 240)
set_x_axis()

plt.sca(production_ax)
plot_single_points(x_superpro[-1], [48], colors.blue_shade.RGBn)
plot_single_points([0, x_superpro[-1]], [141, 70], colors.orange_shade.RGBn)
ethanol_production = get_metric('Ethanol production') / (1e6*rho_etoh)
plot_montecarlo(ethanol_production,
                colors.orange_tint.RGBn,
                colors.orange_shade.RGBn,
                0)
plt.ylim(0, 200)
set_x_axis(False)
plot_vertial_line(1)

plt.sca(production_cost_ax)
plot_single_points(x_superpro[1:], [0.89, 0.84, 0.76], colors.blue_shade.RGBn)
plot_single_points(x_superpro, [0.48, 0.46, 0.44, 0.4], colors.orange_shade.RGBn)
ethanol_production_cost = get_metric('Ethanol production cost') / ethanol_production / 1e6
plot_montecarlo(ethanol_production_cost,
                colors.orange_tint.RGBn,
                colors.orange_shade.RGBn,
                0)
plot_vertial_line(1)
plt.ylim(0, 1.2)
set_x_axis()

plt.subplots_adjust(hspace=.0)
TCI_ax.set_yticks(np.linspace(120, 220, 6))
IRR_ax.tick_params(axis="x", direction="inout", length=4)
IRR_ax.set_zorder(1e6)
production_ax.tick_params(axis="x", direction="inout", length=4)
production_ax.set_zorder(1e6)
production_cost_ax.set_yticks(np.linspace(0, 1, 6))

plt.sca(production_ax)
legend = DoubleColorLegend()
legend.add_box('BioSTEAM')
legend.add_circle('SuperPro (Huang 2016')
legend.legend()