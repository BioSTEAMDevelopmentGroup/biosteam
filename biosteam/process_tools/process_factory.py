# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

class ProcessFactory:
    __slots__ = ('_areas', '_connections', '_inlets', '_outlets', '_names')
    
    def __init__(self):
        self._areas = {}
        self._connections = {}
        self._numbers = {}
        self._inlets = {}
        self._outlets = {}
        self._data = (self._areas, self._connections, self._inlets, self._outlets)
        self._names = set()
    
    def copy(self):
        cls = type(self)
        new = cls.__new__(cls)
        new._areas = self._areas.copy()
        new._connections = self._connections.copy()
        new._numbers = self._numbers.copy()
        new._inlets = self._inlets.copy()
        new._outlets = self._outlets.copy()
        new._data = (new._areas, new._connections, new._numbers, new._inlets, new._outlets)
        new._names = self._names.copy()
        return new
    
    def _add_name(self, name):
        if name in self._names: raise ValueError(f'name {repr(name)} already defined')
        self._names.add(name) 
    
    def add_area(self, name, system_factory, number=None):
        self._add_name(name)
        self._areas[name] = (system_factory, number)
        
    def add_connection(self, name, connection_pipe):
        self._add_name(name)
        self._connections[name] = connection_pipe
        
    def add_inlet(self, name, inlet_pipe):
        self._add_name(name)
        self._inlets[name] = inlet_pipe
    
    def add_outlet(self, name, outlet_pipe):
        self._add_name(name)
        self._outlets[name] = outlet_pipe    
    
    def reset_area(self, name, system_factory, number=None):
        if name not in self._areas: raise ValueError('area {repr(name)} does not exist')
        self._areas[name] = (system_factory, number)
        
    def reset_connection(self, name, connection_pipe):
        if name not in self._connections: raise ValueError('connection {repr(name)} does not exist')
        self._connections[name] = connection_pipe
        
    def reset_inlet(self, name, inlet_pipe):
        if name not in self._inlets: raise ValueError('inlet {repr(name)} does not exist')
        self._inlets[name] = inlet_pipe
    
    def reset_outlet(self, name, outlet_pipe):
        if name not in self._outlets: raise ValueError('outlet {repr(name)} does not exist')
        self._outlets[name] = outlet_pipe    
    
    def __call__(self, ID, **streams):
        areas, connections, numbers, inlets, outlets = self._data
        dct = locals()
        with bst.System(ID) as sys:
            for name, (SF, N) in areas.items(): dct[name] = SF(name, area=N, mockup=True)
            for name, stream in streams.items():
                if name in inlets: eval(inlets[name]).replace(stream)
                elif name in outlets: eval(outlets[name]).replace(stream)
                else: raise ValueError(f'no inlet or outlet named {repr(name)}')
            for connection in connections.values(): eval(connection)
        return sys
    
    def _ipython_display_(self):
        areas, connections, numbers, inlets, outlets = self._data
            
        
        