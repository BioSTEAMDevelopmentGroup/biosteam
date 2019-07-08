# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 22:49:58 2018

@author: yoelr
"""
from .. import Unit, Stream
from .decorators import cost
from ..reaction import Reaction, ParallelReaction

@cost('Volume', 'Reactor', units='m^3',
      CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3,
      fsize=lambda self: self._tau*self._outs[0].volnet/0.8)
class Transesterification(Unit):
    """Create a transesterification reactor that converts 'Lipid' and 'Methanol' to 'Biodiesel' and 'Glycerol'. Finds the amount of catalyst 'NaOCH3' required and consumes it to 'NaOH' and 'Methanol'.
    
    **Parameters**
    
        **efficiency:** Efficiency of conversion (on a 'Lipid' basis)
        
        **methanol2lipid:** Methanol feed to lipid molar ratio
        
        **catalyst_molfrac:** Catalyst to methanol molar ratio
        
        **T:** Operating temperature (K)
    
    **ins**
    
        [0] Lipid feed
        
        [1] Methanol feed (includes catalyst)
        
    **outs**
    
        [0] Product
    
    """
    _kwargs = {'efficiency': None,  # fraction of theoretical conversion
               'methanol2lipid': None,  # Methanol to lipid molar ratio
               'T': None, # operating temperature (K)    
               'catalyst_molfrac': None}  # Catalyst molar fraction in methanol
    _bounds = {'Volume': (0.1, 20)}
    _tau = 1
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1

    def _init(self):
        kw = self._kwargs
        #: [ParallelReaction] Transesterification and catalyst consumption reaction
        self.reaction = ParallelReaction([
          Reaction('Lipid + 3Methanol -> 3Biodiesel + Glycerol',
                   reactant='Lipid',  X=kw['efficiency']),
          Reaction('NaOCH3 -> NaOH + Methanol',
                   reactant='NaOCH3', X=1)])
        self._methanol_composition = Stream.species.kwarray(
                Methanol=1-kw['catalyst_molfrac'],
                NaOCH3=kw['catalyst_molfrac'])
        self._lipid_index, self._methanol_index, self._catalyst_index = \
                Stream.indices(['Lipid', 'Methanol', 'NaOCH3'])
        self._methanol2lipid = kw['methanol2lipid']
        self.T = kw['T'] #: Operation temperature (K).
    
    @property
    def methanol2lipid(self):
        """Methanol feed to lipid molar ratio."""
        return self._methanol2lipid
    @methanol2lipid.setter
    def methanol2lipid(self, ratio):
        self._methanol2lipid = ratio

    @property
    def catalyst_molfrac(self):
        """Catalyst molar fraction in methanol feed."""
        return self._methanol_composition[self._catalyst_index]
    @catalyst_molfrac.setter
    def catalyst_molfrac(self, molfrac):
        meoh = self._methanol_composition
        meoh[self._catalyst_index] = molfrac
        meoh[self._methanol_index] = 1-molfrac

    def _run(self):
        feed, methanol = self.ins
        product, = self.outs
        methanol.mol[:] = (self._methanol_composition
                         * feed.mol[self._lipid_index]
                         * self._methanol2lipid)
        product.mol[:] = feed.mol + methanol.mol
        self.reaction(product.mol)
        product.T = self.T

    def _end(self):
        self._heat_utilities[0](self._Hnet, self.outs[0].T)
        
    
        
        
