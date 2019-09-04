#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 09:36:15 2018

All species for the biodiesel production section of the lipid cane baseline biorefinery are defined here in the Species object, biodiesel_species.

@author: Yoel Rene Cortes-Pena
"""
from biosteam import Species
from biosteam.compounds import Chemical, Substance
from .lipid import lipid

__all__ = ['biodiesel_species']

# %% Make species

sp = Species('Methanol', 'Glycerol', 'Water')
sp.Lipid = lipid # Triolein
sp.Biodiesel = Chemical('methyl oleate')
# sp.FreeLipid = Chemical('oleic acid')

# %% Enthalpy of formations
sp.Lipid.Hf = -2193.7 * 1000 # https://webbook.nist.gov/cgi/cbook.cgi?ID=C122327&Mask=2
#sp.free_lipid.Hfm = -764.80 * 1000
sp.Biodiesel.Hf = -727.64 * 1000
sp.Glycerol.Hf = -669.6 * 1000

# %% Missing properties

sp.Lipid.Tb = 879.9 # Boiling point of Triolein
sp.Lipid.dipole = 0.1
sp.Biodiesel.dipole = 0

# %% Other solids

# Acid and base are basically water
HCl = Substance('HCl', MW=36.46094, rho=10**8)
NaOH = Substance('NaOH', MW=39.997109, rho=10**8)
NaOH.Hfm = 0

# Catalyst is basically methanol
NaOCH3 = Substance('NaOCH3', obj=sp.Methanol, MW=54.023689,
                         rho=1/sp.Methanol.Vm*54.03)
NaOCH3.Hfm = sp.Methanol.Hfm

biodiesel_species = Species.tospecies(tuple(sp) + (HCl, NaOH, NaOCH3))

# Set working species (and in this order)
biodiesel_species.IDs = ('Lipid', 'Methanol', 'Glycerol', 'Biodiesel',
                        'Water', 'NaOH', 'HCl', 'NaOCH3')