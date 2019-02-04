# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:34:07 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam.meta_classes import metaFinal
from biosteam.units.mixer import Mixer


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
    
    """
    kwargs = {'split': None}
    _N_outs = 2

    def __init__(self, ID, outs_ID, ins_ID, **kwargs):
        self.kwargs = kwargs
        self._init_ins(ins_ID)
        self._init_outs(outs_ID)

    def _run(self)
        # Unpack
        split = self.kwargs['split']
        top, bot = self.outs
        if len(self.ins) > 1:
            Mixer.run(self)
        else:
            top.copy_like(self.ins[0])
        bot.copy_like(top)
        top.mol = top.mol*split
        bot.mol = bot.mol - top.mol


class InvSplitter(Unit, metaclass=metaFinal):
    """Create a splitter that sets input stream based on output streams. Must have only one input stream. The output streams are the same temperature, pressure and phase as the input.
    """
    kwargs = {}
    _Graphics = Splitter._Graphics
    
    def __init__(self, ID, outs_ID, ins_ID):
        self._init_ins(ins_ID)
        self._init_outs(outs_ID)
        
    def _run(self)
        feed = self.ins[0]
        feed.mol = self._mol_out
        for out in self.outs:
            out.T, out.P, out.phase = feed.T, feed.P, feed.phase 
