# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import Unit
from ._graphics import hidden_node_graphics
from .utils.piping import HiddenConnection
from .utils import static
from thermosteam import Stream

__all__ = ('hidden_connection',)

def hidden_connection(source, sink):
    upstream = source.outs[0]
    downstream = sink.ins[0]
    stream = Stream(None)
    upstream.sink.ins.replace(upstream, stream)
    hidden_connection = HiddenConnection()
    HiddenConnectionSource(None, ins=upstream, outs=[stream, hidden_connection])
    stream = Stream(None)
    downstream.sink.ins.replace(downstream, stream)
    HiddenConnectionSink(None, ins=[downstream, hidden_connection], outs=stream)
    hidden_connection

@static(N_ins=1, N_outs=2)
class HiddenConnectionSource(Unit):
    _graphics = hidden_node_graphics
    
@static(N_ins=2, N_outs=1) 
class HiddenConnectionSink(Unit):
    _graphics = hidden_node_graphics
