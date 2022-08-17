# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__all__ = ('SystemMesh', 'meshable')

import biosteam as bst
from typing import Sequence

def meshable(obj):
    """Decorator to enable a function or sequence to be accepted to SystemMesh objects."""
    try:
        if not hasattr(obj, 'ins') or not hasattr(obj, 'outs'):
            if isinstance(obj, Sequence):
                obj = bst.MockSystem(obj)
            elif callable(obj):
                if not hasattr(obj, 'ID'): obj.ID = None
                if not hasattr(obj, 'ins'): obj.ins = []
                if not hasattr(obj, 'outs'): obj.outs = []
            else:
                assert False
    except:
        raise TypeError(f"'{obj}' is not meshable; area must be a System, "
                        "MockSystem, SystemFactory, Unit, Sequence[Unit], or function"
                        "object")
    return obj

class SystemMesh:
    __slots__ = ('_areas', '_connections', '_inlets', '_outlets')
    
    def __init__(self):
        self._areas = {}
        self._connections = {}
        self._inlets = {}
        self._outlets = {}
    
    def copy(self):
        cls = type(self)
        new = cls.__new__(cls)
        new._areas = self._areas.copy()
        new._connections = self._connections.copy()
        new._numbers = self._numbers.copy()
        new._inlets = self._inlets.copy()
        new._outlets = self._outlets.copy()
        return new
    
    def add(self, name, obj, number=None, autoconnect=True, **kwargs):
        area = meshable(obj)
        areas = self._areas
        old_area = areas[name] if name in areas else None
        self._areas[name] = (area, number, kwargs)
        if old_area:
            self._inlets.clear()
            self._outlets.clear()
            for i, (j, k, _) in areas.items():
                self._load_area(i, j, k, autoconnect)
        else:
            self._load_area(name, area, number, autoconnect)
        
    def _load_area(self, name, area, number, autoconnect):
        isa = isinstance
        isfunc = callable
        for index, inlet in enumerate(area.ins):
            if isa(inlet, dict):
                ID = inlet.get('ID')
                if ID: self._define_inlet(ID, name, index, autoconnect)
            elif isa(inlet, str):
                ID = inlet
                if ID: self._define_inlet(ID, name, index, autoconnect)
            elif isfunc(inlet):
                ID = inlet.__name__
                self._define_inlet(ID, name, index, autoconnect)
            elif isa(inlet, bst.Stream):
                ID = inlet.ID
                self._define_inlet(ID, name, index, autoconnect)
        for index, outlet in enumerate(area.outs):
            if isa(outlet, dict):
                ID = outlet.get('ID')
                if ID: self._define_outlet(ID, name, index, autoconnect)
            elif isa(outlet, str):
                ID = outlet
                if ID: self._define_outlet(ID, name, index, autoconnect)
            elif isfunc(outlet):
                ID = outlet.__name__
                self._define_outlet(ID, name, index, autoconnect)
            elif isa(outlet, bst.Stream):
                ID = outlet.ID
                self._define_outlet(ID, name, index, autoconnect)
        
    def _rename(self, dct, name):
        other = dct.pop(name)
        try:
            base, num = name.rsplit('_', 1)
            num = int(num)
        except:
            new_name = name + '_1'
        else:
            new_name = base + '_' + str(num + 1) 
        dct[new_name] = other
        base, num = new_name.rsplit('_', 1)
        name = base + '_' + str(int(num) + 1)
        return name
        
    def _define_inlet(self, name, area, index, autoconnect=True):
        if area not in self._areas: raise ValueError('area {repr(area)} does not exist')
        inlets = self._inlets
        value = (str(area), int(index))
        if name in inlets and value != inlets[name]:
            name = self._rename(inlets, name)
        inlets[name] = value
        if autoconnect and name in self._outlets: self.connect(name, name)
    
    def _define_outlet(self, name, area, index, autoconnect=True):
        if area not in self._areas: raise ValueError('area {repr(area)} does not exist')
        outlets = self._outlets
        value = (str(area), int(index))
        if name in outlets and value != outlets[name]:
            name = self._rename(outlets, name)
        outlets[name] = value
        if autoconnect and name in self._inlets: self.connect(name, name)
        
    def connect(self, outlet, inlet):
        if outlet not in self._outlets:
            raise ValueError('outlet {repr(outlet)} does not exist')
        if inlet not in self._inlets:
            raise ValueError('inlet {repr(inlet)} does not exist')
        sink, sink_index = self._inlets[inlet]
        source, source_index = self._outlets[outlet]
        sink = self._areas[sink][0]
        source = self._areas[source][0]
        sink_stream = sink.ins[sink_index]
        source_stream = source.outs[source_index]
        if hasattr(sink_stream, 'sink') and hasattr(source_stream, 'source'):
            sink.ins[sink_index] = source.outs[source_index]
        else:
            self._connections[outlet] = inlet
    
    def __call__(self, ID, **streams):
        areas = self._areas
        connections = self._connections
        inlets = self._inlets 
        outlets = self._outlets
        objs = {}
        with bst.System(ID) as sys:
            for name, (obj, N, kwargs) in areas.items(): 
                if isinstance(obj, bst.SystemFactory):
                    objs[name] = obj(name, area=N, mockup=True, **kwargs)
                elif isinstance(obj, (bst.System, bst.MockSystem)):
                    objs[name] = system = obj
                    if N is not None: # Rename if given area number
                        for i in system.units: i.ID = N
                    # Add units to context management
                    units = system.units
                    bst.main_flowsheet.unit.context_levels[-1].extend(units)
                elif isinstance(obj, bst.Unit):
                    objs[name] = unit = obj
                    if N is not None: unit.ID = N # Rename if given area number
                    # Add units to context management
                    bst.main_flowsheet.unit.context_levels[-1].append(unit)
                elif callable(obj):
                    obj(**kwargs)
                else:
                    raise RuntimeError(f"invalid area type '{type(obj).__name__}'")
            for name, stream in streams.items():
                if name in inlets:
                    objname, index = inlets[name]
                    objs[objname].ins[index] = stream
                elif name in outlets: 
                    objname, index = outlets[name]
                    objs[objname].outs[index] = stream
                else: raise ValueError(f'no inlet or outlet named {repr(name)}')
            for outname, inname in connections.items():
                upstream, index = outlets[outname]
                stream = objs[upstream].outs[index]
                downstream, index = inlets[inname]
                objs[downstream].ins[index] = stream
        for i in sys.units:
            if getattr(i, 'autopopulate', False): i.ins.clear()
        return sys
    
    def show(self):
        areas, connections, inlets, outlets = (self._areas, self._connections, self._inlets, self._outlets)
        info = f"{type(self).__name__}:"
        def get_description(system):
            if isinstance(system, bst.SystemFactory):
                return system.f.__name__
            elif isinstance(system, bst.System):
                return system.ID
            else:
                return None
        fields = {i: {'description': get_description(j), 'number': N, 'ins': {}, 'outs': {}} for i, (j, N, kw) in areas.items()}
        connections = connections.copy()
        for i, j in connections.items(): connections[j] = i
        for name, (area, index) in inlets.items():
            dct = fields[area]['ins']
            if name in connections:
                outname = connections[name]
                outarea, outindex = outlets[outname]
                description = f"{name} from {outarea}-{outindex}"
            else:
                stream = areas[area][0].ins[index]
                if hasattr(stream, 'source') and stream.source:
                    description = f"{name} from {stream.source}-{stream.source.outs.index(stream)}"
                else:
                    description = name
            dct[index] = description
        for name, (area, index) in outlets.items():
            dct = fields[area]['outs']
            if name in connections:
                inname = connections[name]
                inarea, inindex = inlets[inname]
                description = f"{name} to {inindex}-{inarea}"
            else:
                stream = areas[area][0].outs[index]
                if hasattr(stream, 'sink') and stream.sink:
                    description = f"{name} to {stream.sink}-{stream.sink.ins.index(stream)}"
                else:
                    description = name
            dct[index] = description
        keys = ('ins', 'outs')
        section_size = max([len(i) for i in keys]) + 1 # for spaces
        for name, dct in fields.items():
            if dct['number']:
                if dct['description']:
                    info += f"\n{name} ({dct['number']}; {dct['description']})"
                else:
                    info += f"\n{name} ({dct['number']})"
            elif dct['description']:
                info += f"\n{name} ({dct['description']})"
            else:
                info += f"\n{name}"
            for key in keys:
                descriptions = dct[key]
                if descriptions:
                    descriptions = sorted([f"[{i}] {j}" for i, j in descriptions.items()])
                    section = key + (section_size - len(key)) * ' '
                    info += '\n' + section
                    spaces = section_size * ' '
                    info += ('\n' + spaces).join(descriptions)
        print(info)
        
    _ipython_display_ = show
        
        