# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

__all__ = ('Segment',)

class Segment:
    __slots__ = ('sink', 'sink_index', 'source', 'source_index')
    
    def __init__(self, start_stream, end_stream):
        self.sink = sink = start_stream.sink # Segment starts here
        self.sink_index = sink.ins.index(start_stream)
        self.source = source = end_stream.source # Segment ends here
        self.source_index = source.outs.index(end_stream)
        
    def disconnect(self, join_ends=False):
        end_stream = self.source.outs[self.source_index] # Segment includes end stream
        end_unit = end_stream.sink # Segment does not include end unit
        inlet_index = end_unit.ins.index(end_stream)
        end_unit.ins[inlet_index] = None
        if join_ends: 
            # Join the upstream and downstream ends.
            # This maintains the path but without the segment.
            end_unit.ins[inlet_index] = self.sink.ins[self.sink_index]
        self.sink.ins[self.sink_index] = None
            
    def insert(self, stream):
        stream.sink.ins.replace(stream, self.source.outs[self.source_index])
        self.sink.ins[self.sink_index] = stream
        
        
