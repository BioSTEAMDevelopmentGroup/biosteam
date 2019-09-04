# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:11:33 2019

@author: yoelr
"""
__all__ = ['pretreatment_species',
           'biodiesel_species',
           'ethanol_species',
           'lipidcane_species']

from biosteam import Species
from .pretreatment import pretreatment_species
from .biodiesel import biodiesel_species
from .ethanol import ethanol_species

lipidcane_species = Species.tospecies([*pretreatment_species,
                                       *biodiesel_species,
                                       *ethanol_species])