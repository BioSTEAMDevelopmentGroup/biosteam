# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:11:27 2019

@author: yoelr
"""

from biosteam import Species
from biosteam.compounds import Chemical, Gas, Substance

__all__ = ['sugarcane_species']

#%% Set up species

Water = Chemical('Water')        
CO2 = Gas('CO2')
Ethanol = Chemical('Ethanol')
glucose = Substance('Glucose', '492-61-5', MW=180.156, rho=1540.000)
sucrose = Substance('Sucrose', '57-50-1', MW=342.2965, rho=1540.000)
dry_yeast = Substance('DryYeast', 'Yeast', MW=1)
phosphoric_acid = Substance('H3PO4', '7664-38-2', MW=97.994, rho=10**8)
Octane = Chemical('Octane')

#%% Cp Gluclose & sucrose

# Asadi, 2005, Tables, Beet-Sugar Handbook, John Wiley & Sons, Inc (2005), pp.  779-801
# http://onlinelibrary.wiley.com/store/10.1002/9780471790990.oth1/asset/oth1.pdf;jsessionid=1B4B5D3655477F46D578541BE0E4CED2.f03t01?v=1&t=jdewgnrl&s=477661e0dfe5c058191f8f1907d2f63bdaca3e59

# Use heat capacity at T = 60, and saturated with sugar
sucrose.Cp = glucose.Cp = ((2540/4184 - 0.20)*4184 /0.8) * 180.156/1000**2  # J/mol/K
sucrose.Cpm = sucrose.Cp * sucrose.MW
glucose.Cpm = glucose.Cp * glucose.MW
dry_yeast.Cp = dry_yeast.Cpm = 1.2

#%% Cp Phophoric_acid

# https://pubs.acs.org/doi/pdf/10.1021/j150567a017
# HEAT CAPACITY OF PHOSPHORIC ACID SOLUTIONS, 15 TO 80'

# Use heat capacity at T = 59.749 C for 84.81% solution in Water
phosphoric_acid.Cp = (4.184*(0.4624 - 0.1519)/0.8481) * 97.994/1000  # J/mol/K
phosphoric_acid.Cpm = phosphoric_acid.Cp * phosphoric_acid.MW

#%% H_vap (J/kg/k)

glucose.Hvapm = 0
sucrose.Hvapm = 0
phosphoric_acid.Hvapm = 0
dry_yeast.Hvapm = 0

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

# %% Dissolved compounds
dissolved = Species()
for name in ('Ash', 'Cellulose', 'Flocculant', 'Hemicellulose', 'Lignin', 'Solids'):
    setattr(dissolved, name, Substance(name, rho=1540))

dissolved.CaO = Substance('CaO', MW=56.0774, rho=1540)    
dissolved.setprops(dissolved.IDs, 'Hf', 0)
dissolved.setprops(dissolved.IDs, 'T_ref', 298.15)

# %% Heat capacities

# References
# https://www.sciencedirect.com/science/article/pii/0032386182901252
# https://link.springer.com/article/10.1007%2Fs10853-013-7815-6

sugarcane_species = sp = Species.tospecies((CO2, Ethanol, Water,
                                            glucose, sucrose, phosphoric_acid,
                                            Octane, dry_yeast) + tuple(dissolved))
sp.Ash.Cp = sp.Ash.Cpm = 0.09*4184/1000  # from Kumar
sp.CaO.Cp = 1023.88/1000  # from Kumar
sp.CaO.Cpm = 1023.88/1000 * 56.0774
sp.Cellulose.Cp = sp.Cellulose.Cpm = 1364/1000  # from paper
sp.Hemicellulose.Cp = sp.Hemicellulose.Cpm = 1364/1000  # same as cellulose
sp.Flocculant.Cp = sp.Flocculant.Cpm = 4184/1000  # from Kumar
sp.Lignin.Cp = sp.Lignin.Cpm = 1364/1000  # from paper
sp.Solids.Cp = sp.Solids.Cpm = 1100/1000  # common value for solids

sp.setprops(['Ash', 'CaO', 'Cellulose', 'Flocculant',
             'Hemicellulose', 'Lignin', 'Solids'], 'Hvapm', 0)

# %% Heat of combustion

# kJ/kmol
sp.Hemicellulose.Hc = 17000
sp.Cellulose.Hc = 17000
sp.Lignin.Hc = 21000
sp.Glucose.Hc = 2800e3
sp.Sucrose.Hc = 5700e3

sp.IDs = ['Ash', 'Cellulose', 'Hemicellulose', 'Lignin',
          'Glucose', 'Sucrose', 'Solids', 'Water', 'Ethanol',
          'Octane', 'DryYeast', 'H3PO4', 'CaO', 'Flocculant', 'CO2']

# The dry heat values of cellulose, hemicellulose and lignin respectively are 17 MJ/kg (7,320 Btu/lb), 16.63 MJ/kg (7165 Btu/lb) and 21.13 MJ/kg (9,105 Btu/lb) [1]

# [1] Murphy W. K., and K. R. Masters. 1978. Gross heat of combustion of northern red oak (Quercus rubra) chemical components. Wood Sci. 10:139-141.
