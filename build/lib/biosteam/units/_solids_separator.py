# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from .._exceptions import UndefinedCompound
from ._splitter import Splitter, run_split_with_mixing
import biosteam as bst

class SolidsSeparator(Splitter):
    """Create SolidsSeparator object.
    
    **Parameters**
    
        **moisture_content:** Fraction of water in solids
    
    **ins**
    
        [:] Input streams
    
    **outs**
    
        [0] Retentate
        
        [1] Permeate
    
    """
    _N_ins = 2
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 order=None, split, moisture_content):
        Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        #: Moisture content of retentate
        self.mositure_content = moisture_content
        try:
            self._water_index = wi = bst.Stream.species._indexdct['7732-18-5']
        except KeyError:
            raise UndefinedCompound('7732-18-5')
        if self._split[wi] != 0:
            raise ValueError('cannot define water split, only moisture content')
    
    def _run(self):
        run_split_with_mixing(self)
        wi = self._water_index
        retentate, permeate = self.outs
        solids = retentate.massnet
        mc = self.mositure_content
        retentate._mol[wi] = water = (solids * mc/(1-mc))/18.01528
        permeate._mol[wi] -= water
        if permeate._mol[wi] < water:
            raise ValueError(f'not enough water for {repr(self)}')
    
