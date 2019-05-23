# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 16:26:41 2019

@author: yoelr
"""

import lipidcane as lc
import biosteam as bs
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np

solve_IRR = lc.lipidcane_tea.solve_IRR
model = bs.evaluation.Model(ID='IRR',
                            system=lc.lipidcane_sys,
                            metric=solve_IRR)

# %% Simple test

Lipid_cane = lc.Lipid_cane # The feedstock stream
@model.parameter(element=Lipid_cane)
def set_feed_price(feedstock_price): Lipid_cane.price = feedstock_price

model.parameter(lc.set_lipid_fraction, Lipid_cane, kind='coupled')

U34 = bs.find('U34')
@model.parameter(element=U34, kind='coupled')
def set_efficiency(fermentation_efficiency): U34.reset(efficiency=fermentation_efficiency)

samples = np.array([[0.02, 0.85, 0.030],
                    [0.05, 0.90, 0.035],
                    [0.10, 0.95, 0.040]])
model.load_samples(samples)
model.evaluate()
    
# %% Sobol test    
Lipid_cane = lc.Lipid_cane
param = model.parameter(lc.set_lipid_fraction, Lipid_cane, kind='coupled')
bounds = {param.name: [0.02, 0.15]}

original_flow = np.array(Lipid_cane.mol)
@model.parameter(element=Lipid_cane, kind='coupled')
def param(annual_capacity):
    Lipid_cane.mol[:] = original_flow*annual_capacity
bounds[param.name] = [0.8, 1.2]

# Feedstock price
@model.parameter(element=Lipid_cane)
def param(lipid_cane_price):
    Lipid_cane.price = lipid_cane_price
bounds[param.name] = [0.030, 0.040]
    
# Ethanol price
Ethanol = bs.find('Ethanol')
@model.parameter(element=Ethanol)
def param(ethanol_price):
    Ethanol.price = ethanol_price
bounds[param.name] = [0.6, 0.8]

# Biodiesel price
Biodiesel = bs.find('Biodiesel')
@model.parameter(element=Biodiesel)
def param(biodiesel_price):
    Biodiesel.price = biodiesel_price
bounds[param.name] = [1.3, 1.6]

# Crude glycerol price
Crude_glycerol = bs.find('Crude_glycerol')
@model.parameter(element=Crude_glycerol)
def param(crude_glycerol_price):
    Crude_glycerol.price = crude_glycerol_price
bounds[param.name] = [0.2, 0.3]
    
problem = {
'num_vars': 6,
'names':  list(bounds.keys()),
'bounds': list(bounds.values())
}
   
samples = saltelli.sample(problem, 5)    
model.load_samples(samples)
model.evaluate()
Si = sobol.analyze(problem,
                   np.asarray(model.table['IRR']))

    
    
    
    
    
    
    
    
    
    
    