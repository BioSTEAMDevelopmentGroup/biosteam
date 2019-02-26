# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 14:35:52 2019

@author: Guest Group
"""

from lipidcane.pretreatment_species import pretreatment_species
# TODO: Import Clarifier from your file
from biosteam import Stream, Clarifier
import numpy as np

Stream.species = pretreatment_species
feed = Stream('feed', T=372.15,
              Glucose       = 20.1,
              H3PO4         = 1.45,
              Flocculant    = 1.59,
              Lignin        = 501,
              Sucrose       = 120,
              CaO           = 6.99,
              Ash           = 110,
              Cellulose     = 938,
              Hemicellulose = 551,
              Lipid         = 11.4,
              Water         = 2.82e+04)

split = np.array([0.522, 0.522, 0.522, 0.   ,
                  0.   , 0.   , 0.522, 0.   ,
                  0.   , 0.   , 0.   , 0.98 ,
                  0.522])

clarifier = Clarifier('clarifier', ins=feed, split=split)
clarifier.simulate()
print(clarifier.results.table())
