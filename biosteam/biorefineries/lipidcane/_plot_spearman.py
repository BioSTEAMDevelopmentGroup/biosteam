# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 00:00:01 2019

@author: yoelr
"""
import pandas as pd
from biosteam.evaluation.evaluation_tools import plot_spearman

# %%  Plot spearman correlations

# Replacement parameter labels
replacement_labels = {
        'Stream-ethanol price': 'Ethanol price',
        'TEA operating days': 'Operating days',
        'Stream-lipid cane price': 'Lipid cane price',
        'Fermentation-R301 efficiency': 'Fermentation efficiency',
        'Stream-biodiesel price': 'Biodiesel price',
        'Stream-crude glycerol price': 'Crude glycerol price',
        'Stream-lime price': 'Lime price',
        'Stream-lipid cane flow rate': 'Lipid cane flow rate',
        'TEA income tax': 'Income tax',
        'Boiler turbogenerator-BT boiler efficiency': 'Boiler efficiency',
        'Boiler turbogenerator boiler exponent': 'Boiler cost exponent',
        'Boiler turbogenerator boiler base cost': 'Boiler base cost',
        'Boiler turbogenerator turbogenerator electricity rate': 'Turbogenerator electricity rate',
        'Boiler turbogenerator turbogenerator base cost': 'Turbogenerator base cost',
        'Mix tank exponent': 'Mix tank cost exponent',
        'Stream-lipid cane lipid fraction': 'Lipid cane lipid fraction',
        'Transesterification-U61 efficiency': 'Transesterification efficiency',
        'Power utility price': 'Electricity price',
        'Cooling tower base cost': 'Cooling tower base cost'}

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
rhos = pd.read_excel('Spearman correlation lipidcane.xlsx',
                      header=[0], index_col=0).iloc[:, 0]
rhos.name = 'IRR'
# Sort most important
rhos = rhos[(rhos.abs()>0.050)] 

# Plot and fix axis labels
fig, ax = plot_spearman(rhos, top=10)
labels = [item.get_text() for item in ax.get_yticklabels()]
new_labels = [replace_label_text(i) for i in labels]
ax.set_yticklabels(new_labels)