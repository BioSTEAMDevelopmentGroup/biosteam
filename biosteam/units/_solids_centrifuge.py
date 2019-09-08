# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:18:36 2018

@author: yoelr
"""
import numpy as np
from .. import Stream
from .decorators import cost
from ._splitter import Splitter

@cost('Solids loading', cost=68040, CE=567, n=0.50, ub=40, BM=2.03,
      N='Number of centrifuges')
class SolidsCentrifuge(Splitter):
    """Create a solids centrifuge that separates out solids according to user defined split. Assume a continuous scroll solid bowl. 
    
    Parameters
    ----------
    ins
        [:] Input streams
    outs
        [0] Liquid stream
        
        [1] Solids stream
    split: array_like
           Component splits to 0th output stream
    order=None : Iterable[str], defaults to Stream.species.IDs
        Species order of split.
    solids : tuple[str]
             IDs of solids.
    
    References
    ----------
        .. [0] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost Estimation". In Product and Process Design Principles: Synthesis, Analysis, and Evaluation (pp. 481-485). New York: Wiley.
    
    """
    _units = {'Solids loading': 'tonn/hr'}
    _minimum_solids_loading = 2

    def __init__(self, ID='', ins=None, outs=(), *,
                 split, order=None, solids=None):
        super().__init__(ID, ins, outs, split=split, order=order)
        self._solids_index = Stream.indices(solids)
    
    def _design(self):
        mass_solids = 0
        index = self._solids_index
        for s in self.ins:
            mass_solids += s.mass[index]
        ts = np.asarray(mass_solids).sum() # Total solids
        ts *= 0.0011023 # To short tons (2000 lbs/hr)
        self._Design['Solids loading'] = ts
        lb = self._minimum_solids_loading
        if ts < lb:
            self._lb_warning('Solids loading', ts, 'tonn/hr', lb)
    
