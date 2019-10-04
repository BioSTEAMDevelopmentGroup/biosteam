# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 12:59:22 2018

This lipid specie is defined here as tripalmitin.

@author: Yoel Rene Cortes-Pena
"""
import numpy as np
import pandas as pd
from biosteam.compounds import Liquid
import os

__all__ = ['lipid']

MW = 885.432
const_props = {'MW': MW,
               'sigma': 34.7*0.001,  # ChemCad 77F, glyceryl palmitate
               'Tc': 954.1,
               'omega': 1.68620002269745,  # http://www.ethermo.us/Mars1803.htm
               'dipole': 3.14782166481018,  # http://www.ethermo.us/Mars1803.htm
               'Hf': -1844000}  # http://www.ethermo.us/Mars1803.htm
dep_props = {}

# Load material properties
dir_path = os.path.dirname(os.path.realpath(__file__))
data_frame = pd.read_excel(dir_path + '/tripalmitin_liquid.xlsx')
T, rho, mu, nu, Cp, k = np.array(data_frame).transpose()

# Density
density = np.poly1d(np.polyfit(T, rho, 4))
dep_props['rho'] = lambda self: density(self.T)
molar_density = np.poly1d(np.polyfit(T, rho/MW, 4))
dep_props['rhom'] = lambda self: molar_density(self.T)
molar_vol = np.poly1d(np.polyfit(T, MW/rho/1000, 4))
dep_props['Vm'] = lambda self: molar_vol(self.T)

# Dynamic viscosity
dviscosity = np.poly1d(np.polyfit(T, mu, 4))
dep_props['mu'] = lambda self: dviscosity(self.T)

# Kinematic viscosity
kviscosity = np.poly1d(np.polyfit(T, nu, 4))
dep_props['nu'] = lambda self: kviscosity(self.T)

# Heat capacity
molar_heat_capacity = np.poly1d(np.polyfit(T, Cp*MW, 4))
dep_props['Cpm'] = lambda self: molar_heat_capacity(self.T)  # J/mol/K
heat_capacity = np.poly1d(np.polyfit(T, Cp, 4))
dep_props['Cp'] = lambda self: heat_capacity(self.T)  # kJ/kg/K

# Thermal conductivity
conductivity = np.poly1d(np.polyfit(T, k, 4))
dep_props['k'] = lambda self: conductivity(self.T)

# Prandtl number
dep_props['Pr'] = lambda self: self.Cp*self.mu/self.k*1000

# Tripalmitin class
Lipid = Liquid.Factory('Lipid', const_props, dep_props)

# Make the new chemical object
lipid = Lipid('Tripalmitin')
for i, j in const_props.items():
    lipid.__dict__[i] = j
lipid.VaporPressure = lambda T: np.exp(6.596 - 5407.33/T)
lipid.phase_ref = 'l'


