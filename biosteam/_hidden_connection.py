# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from ._graphics import hidden_node_graphics
from .utils.piping import HiddenConnection
from thermosteam import Stream

__all__ = ('hidden_connection',)

def hidden_connection(source, sink):
    upstream = source.outs[0]
    downstream = sink.ins[0]
    hidden_connection = HiddenConnection()
    new_units = []
    if isinstance(upstream.sink, HiddenConnectionSource):
        upstream.sink.outs.append(hidden_connection)
    else:
        stream = Stream(None)
        upstream.sink.ins.replace(upstream, stream)
        new_units.append(
            HiddenConnectionSource(None, ins=upstream, outs=[stream, hidden_connection])
        )
    if isinstance(downstream.source, HiddenConnectionSink):
        downstream.source.ins.append(hidden_connection)
    else:
        stream = Stream(None)
        downstream.sink.ins.replace(downstream, stream)
        new_units.append(
            HiddenConnectionSink(None, ins=[downstream, hidden_connection], outs=stream)
        )
    return new_units

class HiddenConnectionSource(Unit, does_nothing=True):
    _graphics = hidden_node_graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    ticket_name = 'hidden'
    
class HiddenConnectionSink(Unit, does_nothing=True):
    _graphics = hidden_node_graphics
    _ins_size_is_fixed = False
    _N_ins = 2
    _N_outs = 1
    ticket_name = 'hidden'
