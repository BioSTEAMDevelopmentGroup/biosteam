# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:14:01 2018

@author: yoelr
"""
from .. import Unit
from .._exceptions import UndefinedCompound
from .decorators import cost
from .metaclasses import splitter, run_split_with_mixing
import biosteam as bst

@cost('Flow rate', units='kg/hr', cost=1.5e6, CE=541.7,
      n=0.6, S=335e3, kW=2010, BM=2.3)
class CrushingMill(Unit, metaclass=splitter):
    """Create CrushingMill object.
    
    **Parameters**
    
        **moisture_content:** Fraction of water in Baggasse
    
    **ins**
    
        [0] Shredded sugar cane
        
        [1] Recycle water
    
    **outs**
    
        [0] Bagasse
        
        [1] Juice
    
    """
    _N_ins = 2
    _kwargs = {'moisture_content': 0.5}
    def _init(self):
        try:
            self._water_index = wi = bst.Stream.species._indexdct['7732-18-5']
        except KeyError:
            raise UndefinedCompound('7732-18-5')
        if self._split[wi] != 0:
            raise ValueError('cannot define water split, only moisture content')
    
    def _run(self):
        run_split_with_mixing(self)
        wi = self._water_index
        bagasse = self.outs[0]
        juice = self.outs[1]
        solids = bagasse.massnet
        mc = self._kwargs['moisture_content']
        bagasse._mol[wi] = water = (solids * mc/(1-mc))/18.01528
        juice._mol[wi] -= water
        if juice._mol[wi] < water:
            raise ValueError(f'not enough water for {repr(self)}')
    
