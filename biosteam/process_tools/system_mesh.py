# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
                obj.load_inlet_ports([])
                obj.load_outlet_ports([])
            elif callable(obj):
                if not hasattr(obj, 'ID'): obj.ID = obj.__name__
                if not hasattr(obj, 'ins'): obj.ins = []
                if not hasattr(obj, 'outs'): obj.outs = []
            else:
                assert False
    except:
        raise TypeError(f"'{obj}' is not meshable; only System, MockSystem, "
                        "SystemFactory, Unit, Sequence[Unit], or callable "
                        "objects are meshable")
    return obj


class SystemMesh:
    __slots__ = ('_objects', '_connections', '_inlets', '_outlets')
    
    def __init__(self):
        self._objects = {}
        self._connections = {}
        self._inlets = {}
        self._outlets = {}
    
    def copy(self):
        cls = type(self)
        new = cls.__new__(cls)
        new._objects = self._objects.copy()
        new._connections = self._connections.copy()
        new._inlets = self._inlets.copy()
        new._outlets = self._outlets.copy()
        return new
    
    def add(self, name, obj, number=None, autoconnect=True, **kwargs):
        obj = meshable(obj)
        objs = self._objects
        old_obj = objs[name] if name in objs else None
        objs[name] = (obj, number, kwargs)
        if old_obj:
            self._connections.clear()
            self._inlets.clear()
            self._outlets.clear()
            for i, (j, k, _) in objs.items():
                self._load_obj(i, j, k, autoconnect)
        else:
            self._load_obj(name, obj, number, autoconnect)
        
    def remove(self, name, reconnect=True, **kwargs):
        objs = self._objects
        objs.pop(name)
        self._connections.clear()
        self._inlets.clear()
        self._outlets.clear()
        for i, (j, k, _) in objs.items():
            self._load_obj(i, j, k, reconnect)
        
    def _load_obj(self, name, obj, number, autoconnect):
        isa = isinstance
        isfunc = callable
        for index, inlet in enumerate(obj.ins):
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
        for index, outlet in enumerate(obj.outs):
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
        
    def _define_inlet(self, name, obj, index, autoconnect=True):
        if obj not in self._objects: raise ValueError('object with name {repr(obj)} does not exist')
        inlets = self._inlets
        value = (str(obj), int(index))
        if name in inlets and value != inlets[name]:
            name = self._rename(inlets, name)
        inlets[name] = value
        if autoconnect and name in self._outlets: self.connect(name, name)
    
    def _define_outlet(self, name, obj, index, autoconnect=True):
        if obj not in self._objects: raise ValueError('object with name {repr(obj)} does not exist')
        outlets = self._outlets
        value = (str(obj), int(index))
        if name in outlets and value != outlets[name]:
            name = self._rename(outlets, name)
        outlets[name] = value
        if autoconnect and name in self._inlets: self.connect(name, name)
        
    def connect(self, outlet, inlet):
        if outlet not in self._outlets:
            raise ValueError('outlet {repr(outlet)} does not exist')
        if inlet not in self._inlets:
            raise ValueError('inlet {repr(inlet)} does not exist')
        self._connections[outlet] = inlet
    
    def disconnect(self, outlet):
        if outlet not in self._outlets:
            raise ValueError('outlet {repr(outlet)} does not exist')
        self._connections.pop(outlet)
    
    def __call__(self, ID, **streams):
        connections = self._connections
        inlets = self._inlets 
        outlets = self._outlets
        objs = {}
        with bst.System(ID) as sys:
            old_units = []
            for name, (obj, N, kwargs) in self._objects.items(): 
                if isinstance(obj, bst.SystemFactory):
                    objs[name] = obj(name, area=N, mockup=True, **kwargs)
                elif isinstance(obj, (bst.System, bst.MockSystem)):
                    objs[name] = system = obj
                    if N is not None: # Rename if given area number
                        for i in system.units: i.ID = N
                    old_units.extend(system.units)
                    # Add units to context management
                    bst.main_flowsheet.unit.context_levels[-1].extend(system.units)
                elif isinstance(obj, bst.Unit):
                    objs[name] = unit = obj
                    if N is not None: unit.ID = N # Rename if given area number
                    old_units.append(unit)
                    # Add units to context management
                    bst.main_flowsheet.unit.context_levels[-1].append(unit)
                elif callable(obj):
                    obj(**kwargs)
                else:
                    raise RuntimeError(f"invalid object type '{type(obj).__name__}'")
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
            for i in old_units:
                if getattr(i, 'autopopulate', False): i.ins.clear()
        return sys
    
    def show(self):
        objs, connections, inlets, outlets = (self._objects, self._connections, self._inlets, self._outlets)
        info = f"{type(self).__name__}:"
        def get_description(obj):
            if isinstance(obj, bst.SystemFactory):
                return obj.f.__name__
            elif hasattr(obj, 'ID') and isinstance(obj.ID, str):
                return obj.ID
            elif isinstance(obj, bst.MockSystem):
                if len(obj.units) < 4:
                    units = obj.units
                    return ', '.join([i.ID for i in units])
                else:
                    units = obj.units[:4]
                    return ', '.join([i.ID for i in units]) + ', ...'
        fields = {i: {'description': get_description(j), 'number': N, 'ins': {}, 'outs': {}}
                  for i, (j, N, _) in objs.items()}
        connections = connections.copy()
        for i, j in connections.items(): connections[j] = i
        for name, (obj, index) in inlets.items():
            dct = fields[obj]['ins']
            if name in connections:
                outname = connections[name]
                outobj, outindex = outlets[outname]
                description = f"{name} from {outobj}-{outindex}"
            else:
                stream = objs[obj][0].ins[index]
                if hasattr(stream, 'source') and stream.source:
                    description = f"{name} from {stream.source}-{stream.source.outs.index(stream)}"
                else:
                    description = name
            dct[index] = description
        for name, (obj, index) in outlets.items():
            dct = fields[obj]['outs']
            if name in connections:
                inname = connections[name]
                inobj, inindex = inlets[inname]
                description = f"{name} to {inindex}-{inobj}"
            else:
                stream = objs[obj][0].outs[index]
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
        
        