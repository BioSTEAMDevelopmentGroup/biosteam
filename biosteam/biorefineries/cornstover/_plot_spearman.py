# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 00:02:57 2019

@author: yoelr
"""
import pandas as pd
from biosteam import colors
from biosteam.evaluation.evaluation_tools import plot_spearman

# %%  Plot Spearman correlations

# Replacement parameter labels
replacement_labels = {
    'Stream-Ethanol price': 'Ethanol price',
    'TEA operating days': 'Operating days',
    'Stream-cornstover price': 'Cornstover price',
    'Fermentation-R301 efficiency': 'Fermentation efficiency',
    'Stream-cellulase price': 'Cellulase price',
    'Stream-cornstover flow rate': 'Cornstover flow rate',
    'TEA income tax': 'Income tax',
    'Saccharification and co fermentation-R301 saccharification conversion': 'Saccharification conversion',
    'Saccharification and co fermentation-R301 ethanol conversion': 'Fermentation ethanol conversion',
    'Boiler turbogenerator-BT boiler efficiency': 'Boiler efficiency',
    'Boiler turbogenerator boiler base cost': 'Boiler base cost',
    'Boiler turbogenerator turbogenerator base cost': 'Turbogenerator base cost',
    'Pretreatment reactor system base cost': 'Pretreatment reactor base cost',
    'Power utility price': 'Electricity price',
    'Cooling tower base cost': 'Cooling tower base cost',
    'Waste water system cost waste water system base cost': 'Wastewater treatment base cost',
    'Waste water system cost waste water system exponent': 'Wastewater treatment exponent'}


def replace_label_text(label_text):
    """Replace label text for graph."""
    name, distribution = label_text.split(' [')
    lb, mid, ub = eval('[' + distribution)
    if 'efficiency' in name:
        distribution = f" ({lb:.2f}, {mid:.2f}, {ub:.2f})"
    else:
        distribution = f" ({lb:.3g}, {mid:.3g}, {ub:.3g})"
    pos = name.find(' (')
    if pos != -1:
        units = str(name[pos:]).replace('(', '[').replace(')', ']')
        if units == ' [USD/kg]':
            units = ' [$\mathrm{USD} \cdot \mathrm{kg}^{-1}$]'
        elif units == ' [USD/kWhr]':
            units = ' [$\mathrm{USD} \cdot \mathrm{kWhr}^{-1}$]'
        elif units == ' [kg/hr]':
            units = ' [$\mathrm{kg} \cdot \mathrm{hr}^{-1}$]'
        name = name[:pos]
    else:
        units = ''
    if name in replacement_labels:
        name = replacement_labels[name]
    return name + units + distribution

# Get data
rhos = pd.read_excel('Spearman correlation cornstover.xlsx',
                     header=[0], index_col=0).iloc[:, 0]

# Get only important parameters
rhos = rhos[rhos.abs()>0.055] 

# Plot and fix axis labels
fig, ax = plot_spearman(rhos, top=10, name='MESP')
labels = [item.get_text() for item in ax.get_yticklabels()]
new_labels = [replace_label_text(i) for i in labels]
ax.set_yticklabels(new_labels)