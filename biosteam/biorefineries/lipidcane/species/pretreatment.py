#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:26:31 2018

All species for the oil and sugar separation (pretreatment) section of the lipid cane baseline biorefinery are defined here in the Species object, pretreatment_species.

@author: Yoel Rene Cortes-Pena
"""
from biosteam import Species
from biosteam.compounds import Substance
from .ethanol import ethanol_species
from .biodiesel import lipid

__all__ = ['pretreatment_species']

sp = Species()
for name in ('Ash', 'Cellulose', 'Flocculant', 'Hemicellulose', 'Lignin', 'Solids'):
    setattr(sp, name, Substance(name, rho=1540))
    
sp.Lipid = lipid
sp.Lipid.ID = 'Lipid'
sp.CaO = Substance('CaO', MW=56.0774, rho=1540)
sp.setprops(sp.IDs, 'Hf', 0)
sp.setprops(sp.IDs, 'T_ref', 298.15)

# %% Heat capacities

# References
# https://www.sciencedirect.com/science/article/pii/0032386182901252
# https://link.springer.com/article/10.1007%2Fs10853-013-7815-6

pretreatment_species = Species.tospecies(tuple(ethanol_species)+tuple(sp))
sp = pretreatment_species
sp.Ash.Cp = sp.Ash.Cpm = 0.09*4184/1000  # from Kumar
sp.CaO.Cp = 1023.88/1000  # from Kumar
sp.CaO.Cpm = 1023.88/1000 * 56.0774
sp.Cellulose.Cp = sp.Cellulose.Cpm = 1364/1000  # from paper
sp.Hemicellulose.Cp = sp.Hemicellulose.Cpm = 1364/1000  # same as cellulose
sp.Flocculant.Cp = sp.Flocculant.Cpm = 4184/1000  # from Kumar
sp.Lignin.Cp = sp.Lignin.Cpm = 1364/1000  # same as cellulose
sp.Solids.Cp = sp.Solids.Cpm = 1100/1000  # common value for solids

sp.setprops(['Ash', 'CaO', 'Cellulose', 'Flocculant',
             'Hemicellulose', 'Lignin', 'Solids'], 'Hvapm', 0)
sp.IDs =  ('Glucose', 'H3PO4', 'Flocculant',
           'Ethanol', 'Lignin', 'Solids', 'Sucrose',
           'CaO', 'Ash', 'Cellulose', 'Hemicellulose',
           'Lipid', 'Water')

# %% Heat of combustion

# kJ/kmol
sp.Hemicellulose.Hc = 17000
sp.Cellulose.Hc = 17000
sp.Lignin.Hc = 21000
sp.Glucose.Hc = 2800e3
sp.Sucrose.Hc = 5700e3
sp.Lipid.Hc = 36e6


# The dry heat values of cellulose, hemicellulose and lignin respectively are 17 MJ/kg (7,320 Btu/lb), 16.63 MJ/kg (7165 Btu/lb) and 21.13 MJ/kg (9,105 Btu/lb) [1]

# [1] Murphy W. K., and K. R. Masters. 1978. Gross heat of combustion of northern red oak (Quercus rubra) chemical components. Wood Sci. 10:139-141.



