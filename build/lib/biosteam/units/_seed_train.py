# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 23:02:37 2019

@author: yoelr
"""
from .. import Unit

class SeedTrain(Unit):
    
    # Number of stages in series
    N_stages = 5
    
    #: Number of parallel seed trains
    N_trains = 2
    
    #: Cycle time for each batch (hr)
    batch_time = 24
    
    #: Turnaround time for each batch (hr)
    turnarround_time = 12
    
    #: wt % media (e.g. corn steep liquor) in each stage 
    media_loading = 0.50
    
    #: Diammonium phosphate loading in g/L of fermentation broth
    DAP = 0.67 
    
    #: Final fermentor volume (gal)
    maxvol = 2e5
    
    def _run(self): pass

    def _design(self): pass

    def _cost(self): pass

    