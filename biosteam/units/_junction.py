# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 23:29:37 2019

@author: yoelr
"""
from .._unit import Unit
from .metaclasses import final
from .._stream import Stream
from .._mixed_stream import MixedStream
from .._utils import missing_stream, Ins, Outs
from .._flowsheet import find

__all__ = ('Junction',)

def _getstreams(self):
    try:
        upstream, = self._ins
        downstream, = self._outs
    except ValueError as error:
        N_ins = len(self._ins)
        N_outs = len(self._outs)
        if N_ins != 1:
            raise RuntimeError(f'a Junction object must have 1 input stream, not {N_ins}')
        elif N_outs != 1:
            raise RuntimeError(f'a Junction object must have 1 output stream, not {N_outs}')
        else:
            raise error
    return upstream, downstream

# %% Connect between different flowsheets

default_species = lambda upstream, downstream: \
        set(upstream.species._IDs).intersection(downstream.species._IDs)

class Junction(Unit, metaclass=final):
    """Create a Junction object that copies specifications from `upstream`
    to `downstream`. This serves to connect streams with different
    Species object.
    
    **Parameters**
        
        **downstream:** [Stream or str] Flow rate, T, P, and phase information
        will be copied from `upstream` to this stream.
        If None, stream will be missing.
        
        **upstream:** [Stream or str] Stream that will be copied to
        `downstream`. If None, stream will be missing.
        
        **species:** list[str] IDs of species to be passed down.
        If None, all species in common will be passed.
    
    **Examples**
    
        :doc:`Junction Example`
    
    """
    _has_cost = False
    _N_ins = _N_outs = 1
    def __init__(self, upstream=None, downstream='', species=None):
        # Init upstream
        if upstream is None:
            self._ins = Ins(self, (missing_stream,))
        elif isinstance(upstream, Stream):
            self._ins = Ins(self, (upstream,))
        elif isinstance(upstream, str):
            self._ins = Ins(self, (Stream(upstream),))
        
        # Init downstream
        if downstream is None:
            self._outs = Outs(self, (missing_stream,))
        elif isinstance(downstream, Stream):
            self._outs = Outs(self, (downstream,))
        elif isinstance(downstream, str):
            self._outs = Outs(self, (Stream(downstream),))
        
        find.unit[self.ID] = self
        upstream, downstream = _getstreams(self)
        if upstream and downstream:
            self.species = default_species(upstream, downstream)
        else:
            self._species = None
    
    def _run(self):
        upstream, downstream = _getstreams(self)
        try:
            if isinstance(upstream, MixedStream):
                downstream.enable_phases()
                downstream._mol[:, self._downindex] = upstream._mol[:, self._upindex]
            else:
                downstream.disable_phases(upstream._phase)
                downstream._mol[self._downindex] = upstream._mol[self._upindex]
        except Exception as error:
            species = default_species(upstream, downstream)
            if self._species and all([(i in species) for i in self._species]): raise error
            downstream._mol[:] = 0
            self.species = species
            self._run()
        downstream.T = upstream.T
        downstream.P = upstream.P
    simulate = _run
    
    @property
    def ID(self):
        upstream, downstream = _getstreams(self)
        if (upstream or downstream):
            return f"{self._ins[0]} to {self._outs[0]}"
        else:
            return "missing streams"
    _ID = ID
    
    @property
    def species(self):
        return self._species
    @species.setter
    def species(self, IDs):
        upstream, downstream = _getstreams(self)
        self._species = tuple(IDs)
        self._upindex = upstream.indices(*IDs)
        self._downindex = downstream.indices(*IDs)
    
    def _info(self, T, P, flow, fraction):
        info = super()._info(T, P, flow, fraction)
        return info[:info.index(':')+1] + info[info.index('\n'):]
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self.ID}>"


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
