# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for splitting flows.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.splitting.Splitter
.. autoclass:: biosteam.units.splitting.Splitter
.. autoclass:: biosteam.units.splitting.PhaseSplitter 
.. autoclass:: biosteam.units.splitting.FakeSplitter
.. autoclass:: biosteam.units.splitting.ReversedSplitter

"""
from .. import Unit
from .._graphics import splitter_graphics
from thermosteam import separations

__all__ = ('Splitter', 'PhaseSplitter', 'FakeSplitter', 'MockSplitter',
           'ReversedSplitter')

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
       >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
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

    Splits can also be altered after creating the splitter:
        
    .. code-block:: python
       
       >>> S1.split = 0.5
       >>> S1.isplit.show()
       SplitIndexer:
        Water    0.5
        Ethanol  0.5
        
       >>> S1.isplit['Water'] = 1.0
       >>> S1.isplit.show()
       SplitIndexer:
        Water    1
        Ethanol  0.5
        
       >>> S1.split = [0.9, 0.8]
       >>> S1.isplit.show()
       SplitIndexer:
        Water    0.9
        Ethanol  0.8

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
    @split.setter
    def split(self, values):
        split = self.split
        if split is not values:
            split[:] = values
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._isplit = self.thermo.chemicals.isplit(split, order)
        
    def _run(self):
        self.ins[0].split_to(*self.outs, self.split)


class PhaseSplitter(Unit):
    """
    Create a PhaseSplitter object that splits the feed to outlets by phase.
    
    Parameters
    ----------
    ins : stream
        Feed.
    outs : streams
        Outlets.
        
    Notes
    -----
    Phases allocate to outlets in alphabetical order. For example,
    if the feed.phases is 'gls' (i.e. gas, liquid, and solid), the phases
    of the outlets will be 'g', 'l', and 's'.
        
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=10, Ethanol=10)
    >>> feed.vle(V=0.5, P=101325)
    >>> s1 = bst.Stream('s1')
    >>> s2 = bst.Stream('s2')
    >>> PS = bst.PhaseSplitter('PS', feed, [s1, s2])
    >>> PS.simulate()
    >>> PS.show()
    PhaseSplitter: PS
    ins...
    [0] feed
        phases: ('g', 'l'), T: 353.88 K, P: 101325 Pa
        flow (kmol/hr): (g) Water    3.86
                            Ethanol  6.14
                        (l) Water    6.14
                            Ethanol  3.86
    outs...
    [0] s1
        phase: 'g', T: 353.88 K, P: 101325 Pa
        flow (kmol/hr): Water    3.86
                        Ethanol  6.14
    [1] s2
        phase: 'l', T: 353.88 K, P: 101325 Pa
        flow (kmol/hr): Water    6.14
                        Ethanol  3.86
    
    """
    _N_ins = 1
    _N_outs = 2
    _graphics = splitter_graphics
    
    def _run(self):
        separations.phase_split(*self.ins, self.outs)


class MockSplitter(Unit):
    """
    Create a MockSplitter object that does nothing when simulated.
    """
    _graphics = Splitter._graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    
    def _run(self): pass

MockSplitter.line = 'Splitter'
FakeSplitter = MockSplitter    


class ReversedSplitter(Unit):
    """
    Create a splitter that, when simulated, sets the inlet stream based 
    on outlet streams. Must have only one input stream. The outlet streams will
    have the same temperature, pressure and phase as the inlet.
    
    """
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


def reversed_split(inlet, outlets):
    inlet.mol[:] = sum([i.mol for i in outlets])
    T = inlet.T
    P = inlet.P
    phase = inlet.phase
    for out in outlets:
        out.T = T
        out.P = P
        out.phase = phase 

  
    