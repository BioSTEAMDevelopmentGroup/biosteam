# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from .utils import piping

__all__ = ('temporary_connection', 'TemporaryUnit')

temporary_units_dump = []

def temporary_connection(source, sink):
    upstream = source.outs[0]
    downstream = sink.ins[0]
    temporary_stream = piping.TemporaryStream()
    if isinstance(upstream.sink, TemporarySource):
        upstream.sink.outs.append(temporary_stream)
    else:
        stream = piping.TemporaryStream()
        old_connection = upstream.get_connection()
        sink = upstream.sink
        if sink: upstream.sink.ins.replace(upstream, stream)
        TemporarySource(upstream, [stream, temporary_stream], old_connection)
    if isinstance(downstream.source, TemporarySink):
        downstream.source.ins.append(temporary_stream)
    else:
        stream = piping.TemporaryStream()
        old_connection = downstream.get_connection()
        downstream.sink.ins.replace(downstream, stream)
        TemporarySink([downstream, temporary_stream], stream, old_connection)

class TemporaryUnit:
    __slots__ = ('ID', '_ID', 'ins', 'outs', '_ins', '_outs', 'old_connection')
    auxiliary_units = ()
    def __init__(self, ins, outs, old_connection):
        temporary_units_dump.append(self)
        self.ID = self._ID = 'TU' + str(len(temporary_units_dump))
        self.old_connection = old_connection
        self.ins = self._ins = piping.Inlets(
            self, self._N_ins, ins, None, self._ins_size_is_fixed, 5
        )
        self.outs = self._outs = piping.Outlets(
            self, self._N_outs, outs, None, self._outs_size_is_fixed, 5
        )
        
    neighborhood = Unit.neighborhood
    get_downstream_units = Unit.get_downstream_units
    get_upstream_units = Unit.get_upstream_units
    _add_upstream_neighbors_to_set = Unit._add_upstream_neighbors_to_set
    _add_downstream_neighbors_to_set = Unit._add_downstream_neighbors_to_set
    __repr__ = Unit.__repr__
    __str__ = Unit.__str__


class TemporarySource(TemporaryUnit):
    __slots__ = ()
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = True
    _outs_size_is_fixed = False
    

class TemporarySink(TemporaryUnit):
    __slots__ = ()
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True
    
    