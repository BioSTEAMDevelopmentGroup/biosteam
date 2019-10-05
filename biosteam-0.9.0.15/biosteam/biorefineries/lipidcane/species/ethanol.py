#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 22:35:30 2018

All species for the ethanol production section of the lipid cane baseline biorefinery are defined here in the Species object, ethanol_species.

@author: Yoel Rene Cortes-Pena
"""
from biosteam import Species
from biosteam.compounds import Chemical, Gas, StaticChemical, Substance

__all__ = ['ethanol_species']

#%% Set up species

Water = Chemical('Water')        
CO2 = Gas('CO2')
Ethanol = Chemical('Ethanol')
_glucosoe = Chemical('Glucose')
glucose = StaticChemical('Glucose')
sucrose = StaticChemical('Sucrose')
dry_yeast = Substance('DryYeast', 'Yeast', MW=1)
phosphoric_acid = Substance('H3PO4', '7664-38-2', MW=97.994, rho=10**8)
Octane = Chemical('Octane')
ethanol_species = sp = Species.tospecies((CO2, Ethanol, Water, glucose, sucrose, phosphoric_acid, Octane, dry_yeast))

#%% Cp Gluclose & sucrose

# # Asadi, 2005, Tables, Beet-Sugar Handbook, John Wiley & Sons, Inc (2005), pp.  779-801
# # http://onlinelibrary.wiley.com/store/10.1002/9780471790990.oth1/asset/oth1.pdf;jsessionid=1B4B5D3655477F46D578541BE0E4CED2.f03t01?v=1&t=jdewgnrl&s=477661e0dfe5c058191f8f1907d2f63bdaca3e59

# # Use heat capacity at T = 60, and saturated with sugar
# sucrose.Cp = glucose.Cp = ((2540/4184 - 0.20)*4184 /0.8) * 180.156/1000**2  # J/mol/K
# sucrose.Cpm = sucrose.Cp * sucrose.MW
# glucose.Cpm = glucose.Cp * glucose.MW
dry_yeast.Cp = dry_yeast.Cpm = 1.2

#%% Cp Phophoric_acid

# https://pubs.acs.org/doi/pdf/10.1021/j150567a017
# HEAT CAPACITY OF PHOSPHORIC ACID SOLUTIONS, 15 TO 80'

# Use heat capacity at T = 59.749 C for 84.81% solution in Water
phosphoric_acid.Cp = (4.184*(0.4624 - 0.1519)/0.8481) * 97.994/1000  # J/mol/K
phosphoric_acid.Cpm = phosphoric_acid.Cp * phosphoric_acid.MW

#%% Enthalpy of formation (all information from NIST)

# The heats of formation (kJ/kmol)
Ethanol.Hf = -277.690 * 1000
Water.Hf = -241.820 * 1000
CO2.Hf = -393.520 * 1000  # in gas
glucose.Hf = -1274*1000
sucrose.Hf = -2221.2*1000
Octane.Hf = -250*1000
phosphoric_acid.Hf = -1271.66*1000
dry_yeast.Hf = 0
