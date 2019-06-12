# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 14:34:07 2018

@author: yoelr
"""
from .. import Unit, Stream
from .metaclasses import final
from numpy import asarray
from .metaclasses._splitter import split

class Splitter(Unit, metaclass=final):
    """Create a splitter that separates mixed streams based on splits.

    **Parameters**

        **split:** Shoulds be one of the following:
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
    results = None
    _has_cost = False
    _N_outs = 2

    def __init__(self, ID='', outs=(), ins=None, split=None, order=None):
        self.ID = ID
        try:
            self._reorder_ = Stream._cls_species._reorder
        except AttributeError:
            raise RuntimeError('must specify Stream.species first')
        self._split = self._reorder_(split, order) if order else asarray(split)
        self._init_ins(ins)
        self._init_outs(outs)

    def _run(self):
        top, bot = self._outs
        feed = self._ins[0]
        net_mol = feed.mol
        bot.T = top.T = feed.T
        bot.P = top.P = feed.P
        bot._phase = top._phase = feed._phase
        top._mol[:] = net_mol * self._split
        bot._mol[:] = net_mol - top._mol
    
    split = split
    simulate = _run
    summary = Unit._cost


class InvSplitter(Unit, metaclass=final):
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
