# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:34:07 2018

@author: yoelr
"""
from .. import Unit, Stream
from .._meta_final import metaFinal

class Splitter(Unit, metaclass=metaFinal):
    """Create a splitter that separates mixed streams based on splits.

    **Parameters**

        **split:** Shoulds be one of the following:
            * *[float]* The fraction of net feed in the 0th output stream
            * *[array_like]* Component wise split of feed to 0th output stream
    
    **ins**
    
        [:] Input streams
        
    **outs**
    
        [0] Split stream
        
        [1] Remainder stream
    
    **Examples**
    
        Create a Splitter object with an ID, any number of input streams, two output streams, and an overall split:
            
        .. code-block:: python
        
           >>> from biosteam import Species, Stream, Splitter
           >>> Stream.species = Species('Water', 'Ethanol')
           >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
           >>> S1 = Splitter('S1', ins=feed, outs=('out1', 'out2'), split=0.1)
           >>> S1.simulate()
           >>> S1.show()
           Splitter: S1
           ins...
           [0] feed
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    20
                               Ethanol  10
           outs...
           [0] out1
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    2
                               Ethanol  1
           [1] out2
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    18
                               Ethanol  9
          
        Create a Splitter object, but this time with a component wise split:
            
        .. code-block:: python
        
           >>> from biosteam import Species, Stream, Splitter
           >>> Stream.species = Species('Water', 'Ethanol')
           >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
           >>> S1 = Splitter('S1', ins=feed, outs=('out1', 'out2'),
           ...               split=(0.1, 0.99))
           >>> S1.simulate()
           >>> S1.show()
           Splitter: S1
           ins...
           [0] feed
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    20
                               Ethanol  10
           outs...
           [0] out1
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    2
                               Ethanol  9.9
           [1] out2
               phase: 'l', T: 340 K, P: 101325 Pa
               flow (kmol/hr): Water    18
                               Ethanol  0.1
    
    """
    results = None
    _has_cost = False
    _N_outs = 2

    def __init__(self, ID='', outs=(), ins=None, split=None, order=None):
        self.ID = ID
        self._reorder = Stream._cls_species._reorder
        self._kwargs = {'split': self._reorder(split, order) if order else split}
        self._init_ins(ins)
        self._init_outs(outs)

    def _reset(self, split, order=None):
        self._kwargs['split'] = self._reorder(split, order) if order else split

    def _run(self):
        split = self._kwargs['split']
        top, bot = self._outs
        ins = self._ins
        if len(ins) > 1: Stream.sum(top, ins)
        else: top.copylike(ins[0])
        bot.copylike(top)
        top._mol[:]*= split
        bot._mol[:]-= top._mol
    
    simulate = _run
    summary = Unit._cost


class InvSplitter(Unit, metaclass=metaFinal):
    """Create a splitter that sets the input stream based on output streams. Must have only one input stream. The output streams will become the same temperature, pressure and phase as the input.
    """
    line = 'Splitter'
    results = None
    _graphics = Splitter._graphics
    _has_cost = False
    
    def __init__(self, ID='', outs=(), ins=None):
        self.ID = ID
        self._init_ins(ins)
        self._init_outs(outs)
        
    def _run(self):
        feed = self.ins[0]
        feed._mol[:] = self._mol_out
        T = feed.T
        P = feed.P
        phase = feed.phase
        for out in self._outs:
            out.T = T
            out.P = P
            out.phase = phase 
            
    simulate = _run
    summary = Unit._cost
