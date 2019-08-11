# -*- coding: utf-8 -*-
"""
Created on Mon May 20 22:04:02 2019

@author: yoelr
"""
from .. import Unit, Stream
from numpy import asarray
from pandas import Series

__all__ = ('run_split', 'run_split_with_mixing')


def run_split(self):
    """Splitter mass and energy balance function based on one input stream."""
    top, bot = self.outs
    feed = self._ins[0]
    net_mol = feed.mol
    top._mol[:] = net_mol * self._split
    bot._mol[:] = net_mol - top._mol
    bot.T = top.T = feed.T
    bot.P = top.P = feed.P

def run_split_with_mixing(self):
    """Splitter mass and energy balance function with mixing all input streams."""
    top, bot = self._outs
    ins = self._ins
    if len(ins) > 1: Stream.sum(top, ins)
    else: top.copylike(ins[0])
    bot.copylike(top)
    top._mol[:] *= self._split
    bot._mol[:] -= top._mol


class Splitter(Unit):
    """Create a splitter that separates mixed streams based on splits.

    **Parameters**

        **split:** Should be one of the following
            * [float] The fraction of net feed in the 0th output stream
            * [array_like] Componentwise split of feed to 0th output stream
            
        **order:** Iterable[str] Species order of split
    
    **ins**
    
        [0] Feed stream
        
    **outs**
    
        [0] Split stream
        
        [1] Remainder stream
    
    **Examples**
    
        Create a Splitter object with an ID, a feed stream, two output streams, and an overall split:
            
        .. code-block:: python
        
           >>> from biosteam import Species, Stream, units
           >>> Stream.species = Species('Water', 'Ethanol')
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
          
        Create a Splitter object, but this time with a componentwise split:
            
        .. code-block:: python
        
           >>> from biosteam import Species, Stream, units
           >>> Stream.species = Species('Water', 'Ethanol')
           >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
           >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
           ...                     split=(0.1, 0.99))
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
        
           >>> from biosteam import Species, Stream, units
           >>> Stream.species = Species('Water', 'Ethanol')
           >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
           >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
           ...                     split=(0.99, 0.01), order=('Ethanol', 'Water'))
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
    
    @property
    def split(self):
        """[array] Componentwise split of feed to 0th output stream."""
        return self._split_series
    
    def __init__(self, ID='', ins=None, outs=(), *, split, order=None):
        Unit.__init__(self, ID, ins, outs)
        self._split = split = Stream.species.array(order, split) if order else asarray(split)
        self._split_series =  Series(split, Stream.species.IDs)

    def _run(self):
        top, bot = self._outs
        feed = self._ins[0]
        net_mol = feed.mol
        bot.T = top.T = feed.T
        bot.P = top.P = feed.P
        bot._phase = top._phase = feed._phase
        top._mol[:] = net_mol * self._split
        bot._mol[:] = net_mol - top._mol


class InvSplitter(Unit):
    """Create a splitter that sets the input stream based on output streams. Must have only one input stream. The output streams will become the same temperature, pressure and phase as the input.
    """
    line = 'Splitter'
    _graphics = Splitter._graphics
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
    


        
    