# -*- coding: utf-8 -*-
"""
Created on Mon May 20 22:04:02 2019

@author: yoelr
"""
from .. import Unit
from .._graphics import splitter_graphics
import numpy as np 
from thermosteam.indexer import ChemicalIndexer
from ._process_specification import ProcessSpecification

__all__ = ('run_split', 'run_split_with_mixing', 'Splitter', 'FakeSplitter',
           'ReversedSplitter', 'ReversedSplit')


def run_split(self):
    """Splitter mass and energy balance function based on one input stream."""
    top, bot = self.outs
    feed = self.ins[0]
    top_mol = top.mol
    feed_mol = feed.mol
    top_mol[:] = feed_mol * self.split
    bot.mol[:] = feed_mol - top_mol
    bot.T = top.T = feed.T
    bot.P = top.P = feed.P

def run_split_with_mixing(self):
    """Splitter mass and energy balance function with mixing all input streams."""
    top, bot = self.outs
    ins = self.ins
    top.mix_from(ins)
    bot.copy_like(top)
    top_mol = top.mol
    top_mol[:] *= self.split
    bot.mol[:] -= top_mol


class Splitter(Unit):
    """
    Create a splitter that separates mixed streams based on splits.

    Parameters
    ----------
    ins : stream
        Inlet fluid to be split.
    outs : stream sequence
        * [0] Split stream
        * [1] Remainder stream    
    split : Should be one of the following
            * [float] The fraction of net feed in the 0th outlet stream
            * [array_like] Componentwise split of feed to 0th outlet stream
            * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    
    Examples
    --------
    Create a Splitter object with an ID, a feed stream, two outlet streams,
    and an overall split:
        
    .. code-block:: python
    
       >>> from biosteam import units, settings, Stream
       >>> settings.set_thermo(['Water', 'Ethanol'])
       >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
       >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'), split=0.1)
       >>> S1.simulate()
       >>> S1.show()
       Splitter: S1
       ins...
       [0] feed
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    20
                           Ethanol  10
       outs...
       [0] top
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    2
                           Ethanol  1
       [1] bot
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    18
                           Ethanol  9
      
    Create a Splitter object, but this time with a componentwise split
    using a dictionary:
        
    .. code-block:: python
    
       >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
       ...                     split={'Water': 0.1, 'Ethanol': 0.99})
       >>> S1.simulate()
       >>> S1.show()
       Splitter: S1
       ins...
       [0] feed
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    20
                           Ethanol  10
       outs...
       [0] top
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    2
                           Ethanol  9.9
       [1] bot
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    18
                           Ethanol  0.1
                           
    Create a Splitter object using componentwise split, but this time specify the order:
        
    .. code-block:: python
    
       >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
       ...                     order=('Ethanol', 'Water'),
       ...                     split=(0.99, 0.10))
       >>> S1.simulate()
       >>> S1.show()
       Splitter: S1
       ins...
       [0] feed
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    20
                           Ethanol  10
       outs...
       [0] top
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    2
                           Ethanol  9.9
       [1] bot
           phase: 'l', T: 340 K, P: 101325 Pa
           flow (kmol/hr): Water    18
                           Ethanol  0.1

    """
    _N_outs = 2
    _graphics = splitter_graphics
    
    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit._data
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        chemicals = self.thermo.chemicals
        if isinstance(split, dict):
            assert not order, "cannot pass 'other' key word argument when split is a dictionary"
            order, split = zip(*split.items())
        
        if order:
            isplit = chemicals.iarray(order, split)
        elif hasattr(split, '__len__'):
            isplit = ChemicalIndexer.from_data(np.asarray(split),
                                               phase=None,
                                               chemicals=chemicals)
        else:
            split = split * np.ones(chemicals.size)
            isplit = ChemicalIndexer.from_data(split,
                                               phase=None,
                                               chemicals=chemicals)
        self._isplit = isplit

    def _run(self):
        top, bot = self.outs
        feed, = self.ins
        feed_mol = feed.mol
        bot.T = top.T = feed.T
        bot.P = top.P = feed.P
        bot.phase = top.phase = feed.phase
        top.mol[:] = top_mol = feed_mol * self.split
        bot.mol[:] = feed_mol - top_mol


class FakeSplitter(Unit):
    """
    Create a FakeSplitter object that does nothing when simulated.
    """
    _graphics = Splitter._graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    
    def _run(self): pass
    
    def create_reversed_splitter_process_specification(self, ID='', ins=None, outs=(),
                                                       description=None):
        return ProcessSpecification(ID, ins, outs, self.thermo, 
                                    specification=ReversedSplit(self),
                                    description=description or 'Reverse split')    
    
FakeSplitter.line = 'Splitter'

class ReversedSplitter(Unit):
    """Create a splitter that, when simulated, sets the inlet stream based on outlet streams. Must have only one input stream. The outlet streams will become the same temperature, pressure and phase as the input."""
    _graphics = Splitter._graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    power_utility = None
    heat_utilities = ()
    results = None
    
    def _run(self):
        inlet, = self.ins
        outlets = self.outs
        reversed_split(inlet, outlets)

class ReversedSplit:
    """Create a splitter that, when called, sets the inlet stream based on outlet streams. Must have only one input stream. The outlet streams will become the same temperature, pressure and phase as the input."""
    __slots__ = ('splitter',)
    
    def __init__(self, splitter):
        self.splitter = splitter
    
    @property
    def __name__(self):
        return 'ReversedSplit'
    
    def __call__(self):
        splitter = self.splitter
        inlet, = splitter.ins
        outlets = splitter.outs
        reversed_split(inlet, outlets)
        
def reversed_split(inlet, outlets):
    inlet.mol[:] = sum([i.mol for i in outlets])
    T = inlet.T
    P = inlet.P
    phase = inlet.phase
    for out in outlets:
        out.T = T
        out.P = P
        out.phase = phase 

  
    