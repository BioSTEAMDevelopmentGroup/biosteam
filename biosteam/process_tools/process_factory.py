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
    
    def _add_name(self, name, safe):
        if safe and name in self._names: raise ValueError(f'name {repr(name)} already defined')
        self._names.add(name) 
    
    def add_area(self, name, system_factory, number=None, safe=True):
        self._add_name(name, safe)
        self._areas[name] = (system_factory, number)
        
    def add_connection(self, name, outlet, inlet, safe=True):
        self._add_name(name, safe)
        outlet = self._outlets.pop(outlet)
        inlet = self._inlets.pop(inlet)
        self._connections[name] = (outlet, inlet)
        
    def add_inlet(self, name, inlet_pipe, safe=True):
        self._add_name(name, safe)
        area, index = inlet_pipe.split('-')
        if name not in self._areas: raise ValueError('area {repr(name)} does not exist')
        self._inlets[name] = (area, int(index))
    
    def add_outlet(self, name, outlet_pipe, safe=True):
        self._add_name(name, safe)
        area, index = outlet_pipe.split('-')
        if name not in self._outlets: raise ValueError('outlet {repr(name)} does not exist')
        self._outlets[name] = (area, int(index))
    
    def __call__(self, ID, **streams):
        areas, connections, numbers, inlets, outlets = self._data
        systems = {}
        with bst.System(ID) as sys:
            for name, (SF, N) in areas.items(): systems[name] = SF(name, area=N, mockup=True)
            for name, stream in streams.items():
                if name in inlets:
                    sysname, index = inlets[name]
                    systems[sysname].ins[index] = stream
                elif name in outlets: 
                    sysname, index = outlets[name]
                    systems[sysname].outs[index] = stream
                else: raise ValueError(f'no inlet or outlet named {repr(name)}')
            for connection in connections.values():
                upstream, uindex = outlets[name]
                downstream, dindex = inlets[name]
                systems[downstream].ins[dindex] = systems[upstream].outs[uindex]
        return sys
    
    def _ipython_display_(self):
        areas, connections, numbers, inlets, outlets = self._data
            
        
        