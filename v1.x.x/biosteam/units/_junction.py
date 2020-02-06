# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 23:29:37 2019

@author: yoelr
"""
from .._unit import Unit
from .._stream import Stream
from .._mixed_stream import MixedStream
from ..utils import MissingStream, Ins, Outs

__all__ = ('Junction',)

def _getstreams(self):
    try:
        upstream, = self._ins
        downstream, = self._outs
    except ValueError as error:
        N_ins = len(self._ins)
        N_outs = len(self._outs)
        print(self._ins)
        print(self._outs)
        if N_ins != 1:
            raise RuntimeError(f'a Junction object must have 1 input stream, not {N_ins}')
        elif N_outs != 1:
            raise RuntimeError(f'a Junction object must have 1 output stream, not {N_outs}')
        else:
            raise error
    return upstream, downstream

# %% Connect between different flowsheets

default_species = lambda upstream, downstream: \
    [i.ID for i in set(upstream.species._compounds).intersection(downstream.species._compounds)]

class Junction(Unit):
    """Create a Junction object that copies specifications from `upstream`
    to `downstream`. This serves to connect streams with different
    Species object.
    
    Parameters
    ----------    
    upstream=None : Stream or str, defaults to missing stream
        Stream that will be copied to `downstream`.
    downstream="" : Stream or str, defaults to new stream
        Flow rate, T, P, and phase information
        will be copied from `upstream` to this stream.
        If None, stream will be missing.
    species=None : list[str], defaults to all species in common
        IDs of species to be passed down.
    
    Examples
    --------
    Create a Junction object and connect streams with different Species objects:
        
    >>> from biosteam import *
    >>> species = Species('Water')
    >>> species.Sugar = compounds.Substance('Sugar')
    >>> Stream.species = species
    >>> s1 = Stream('s1', Water=20)
    >>> Stream.species = Species('Water', 'Ethanol')
    >>> s2 = Stream('s2') # Note that s2 and s1 have different Species objects
    >>> J1 = units.Junction(upstream=s1, downstream=s2)
    >>> J1.simulate()
    >>> J1.show()
    Junction:
    ins...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  20
    outs...
    [0] s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  20
        
    """
    _has_cost = False
    _N_ins = _N_outs = 1
    def __init__(self, ID="", upstream=None, downstream='', species=None):
        # Init upstream
        if upstream is None:
            self._ins = Ins(self, (MissingStream,))
        elif isinstance(upstream, Stream):
            self._ins = Ins(self, (upstream,))
        elif isinstance(upstream, str):
            self._ins = Ins(self, (Stream(upstream),))
        
        # Init downstream
        if downstream is None:
            self._outs = Outs(self, (MissingStream,))
        elif isinstance(downstream, Stream):
            self._outs = Outs(self, (downstream,))
        elif isinstance(downstream, str):
            self._outs = Outs(self, (Stream(downstream),))
        
        self.species = species
        self.ID = ID
    
    def _material_balance(self, upstream, downstream):
        if isinstance(upstream, MixedStream):
            downstream.enable_phases()
            downstream._mol[:, self._downindex] = upstream._mol[:, self._upindex]
        else:
            downstream.disable_phases(upstream._phase)
            downstream._mol[self._downindex] = upstream._mol[self._upindex]
    
    def _run(self):
        upstream, downstream = _getstreams(self)
        if not self._species:
            self.species = default_species(upstream, downstream)
        self._material_balance(upstream, downstream)
        downstream.T = upstream.T
        downstream.P = upstream.P
    simulate = _run
    
    @property
    def species(self):
        return self._species
    @species.setter
    def species(self, IDs):
        if IDs is None:
            self._species = IDs
        else:
            upstream, downstream = _getstreams(self)
            self._species = IDs = tuple(IDs)
            self._upindex = upstream.indices(IDs)
            self._downindex = downstream.indices(IDs)


def node_function(self):
    node = self._graphics.node
    if not any(_getstreams(self)):
        node['fontsize'] = '18'
        node['shape'] = 'plaintext'
        node['fillcolor'] = 'none'
    else:
        node['width'] = '0.1'
        node['shape'] = 'point'
        node['color'] = node['fillcolor'] = 'black'

Junction._graphics.node_function = node_function
del node_function
