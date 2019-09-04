# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from biosteam import colors
from biosteam.utils import DoubleColorLegend
from biosteam.evaluation.evaluation_tools import plot_spearman, plot_single_points, \
     plot_montecarlo_across_coordinate, plot_montecarlo

    
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

# %%  Plot Spearman correlations

# Replacement parameter labels
replacement_labels = {
    'Stream-Ethanol price': 'Ethanol price',
    'TEA operating days': 'Operating days',
    'Stream-Cornstover price': 'Cornstover price',
    'Fermentation-U34 efficiency': 'Fermentation efficiency',
    'Stream-Cellulase422 price': 'Cellulase price',
    'Stream-Cornstover flow rate': 'Cornstover flow rate',
    'TEA income tax': 'Income tax',
    'Saccharification and co fermentation-U22 saccharification conversion': 'Saccharification conversion',
    'Saccharification and co fermentation-U22 ethanol conversion': 'Ethanol conversion',
    'Boiler turbogenerator-BT boiler efficiency': 'Boiler efficiency',
    'Boiler turbogenerator boiler base cost': 'Boiler base cost',
    'Boiler turbogenerator turbogenerator base cost': 'Turbogenerator base cost',
    'Pretreatment reactor system base cost': 'Pretreatment reactor base cost',
    'Power utility price': 'Electricity price',
    'Cooling tower base cost': 'Cooling tower base cost',
    'Waste water system cost waste water system base cost': 'Wastewater treatment base cost'}


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
rhos = pd.read_excel('Spearman correlation cornstover.xlsx',
                     header=[0], dtype=float)['Minimum ethanol selling price']

# Get only important parameters
rhos = rhos[rhos.abs()>0.055] 

# Plot and fix axis labels
fig, ax = plot_spearman(rhos, top=10)
labels = [item.get_text() for item in ax.get_yticklabels()]
new_labels = [replace_label_text(i) for i in labels]
ax.set_yticklabels(new_labels)
