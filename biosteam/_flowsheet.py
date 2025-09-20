# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
As BioSTEAM objects are created, they are automatically registered. 
The `main_flowsheet` object allows the user to find any Unit, Stream or System instance.
When `main_flowsheet` is called, it simply looks up the item and returns it. 
"""
from __future__ import annotations
from typing import Optional, Iterable
import biosteam as bst
from thermosteam.utils import Registry
from thermosteam import AbstractStream
from ._unit import AbstractUnit
from ._system import System

__all__ = ('main_flowsheet', 'Flowsheet', 'F')

# %% Flowsheet search      

class FlowsheetRegistry:
    __getitem__ = object.__getattribute__
    
    def clear(self):
        self.__dict__.clear()
        main_flowsheet.set_flowsheet('default')
    
    def __contains__(self, obj):
        dct = self.__dict__
        if isinstance(obj, str):
            return obj in dct
        elif isinstance(obj, Flowsheet):
            return obj.ID in dct and dct[obj.ID] is obj
        else:
            return False
    
    def __setattr__(self, key, value):
        raise TypeError(f"'{type(self).__name__}' object does not support attribute assignment")
    __setitem__ = __setattr__
    
    def __iter__(self):
        return self.__dict__.values().__iter__()
    
    def __delattr__(self, key):
        if key == main_flowsheet.ID:
            raise AttributeError('cannot delete main flowsheet')
        else:
            super().__delattr__(key)
            
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join([str(i) for i in self])}>"
    
    def _ipython_display_(self): # pragma: no cover
        print(f'{type(self).__name__}:\n ' + '\n '.join([str(i) for i in self]))
    
    
class Flowsheet:
    """
    Create a Flowsheet object which stores references to all stream, unit,
    and system objects. For a tutorial on flowsheets, visit
    :doc:`../tutorial/Managing_flowsheets`.
	
    """
    
    line: str = "Flowsheet"
    
    #: All flowsheets.
    flowsheet: FlowsheetRegistry = FlowsheetRegistry()
    
    def __new__(cls, ID=None):        
        self = super().__new__(cls)
        if ID is None: ID = ''
        
        #: Contains all System objects as attributes.
        self.system: Registry = Registry()
        
        #: Contains all Unit objects as attributes.
        self.unit: Registry = Registry()
        
        #: Contains all Stream objects as attributes.
        self.stream: Registry = Registry()
        
        #: ID of flowsheet.
        self._ID: str = ID
        
        #: Temporary flowsheet stack.
        self._temporary_stack = []
        
        self.flowsheet.__dict__[ID] = self
        return self
    
    def __enter__(self):
        """
        Temporarily register all objects in this flowsheet instead of 
        the main flowsheet.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> f = bst.Flowsheet('f')
        >>> with f:
        ...     M1 = bst.Mixer('M1')
        >>> M1 in bst.main_flowsheet.unit
        False
        >>> M1 in f.unit
        True
        
        """
        self._temporary_stack.append(main_flowsheet.get_flowsheet())
        main_flowsheet.set_flowsheet(self)
        return self
    
    def __exit__(self, type, exception, traceback):
        main_flowsheet.set_flowsheet(self._temporary_stack.pop())
        if exception: raise exception
    
    def temporary(self): return self # for backwards compatibility
    
    def __reduce__(self):
        return self.from_registries, self.registries
    
    def __getattr__(self, name):
        obj = (self.stream.search(name)
               or self.unit.search(name)
               or self.system.search(name))
        if not obj: raise AttributeError(f"no registered item '{name}'")
        return obj
    
    def __setattr__(self, key, value):
        if self in self.flowsheet.__dict__:
            raise AttributeError("cannot register object through flowsheet")
        else:
            super().__setattr__(key, value)
    
    @property
    def ID(self):
        return self._ID
    
    @classmethod
    def from_registries(cls, ID, stream, unit, system):
        flowsheet = super().__new__(cls)
        flowsheet.stream = stream
        flowsheet.unit = unit
        flowsheet.system = system
        flowsheet._ID = ID
        flowsheet._temporary_stack = []
        flowsheet.flowsheet.__dict__[ID] = flowsheet
        return flowsheet
    
    @property
    def registries(self):
        return (self.stream, self.unit, self.system)
    
    def clear(self, reset_ticket_numbers=True):
        for registry in self.registries: registry.clear()
        if reset_ticket_numbers:
            for i in (AbstractStream, AbstractUnit, System): i.ticket_numbers.clear()
    
    def discard(self, ID):
        for registry in self.registries: registry.discard(ID)
    
    def remove_unit_and_associated_streams(self, ID):
        stream_registry = self.stream
        unit = self.unit.pop(ID)
        for inlet in unit._ins:
            if inlet.source: continue
            stream_registry.discard(inlet)
        for outlet in unit._outs:
            if outlet._sink: continue
            stream_registry.discard(outlet)
    
    def update(self, flowsheet):
        for registry, other_registry in zip(self.registries, flowsheet.registries):
            registry.data.update(other_registry.data)
    
    def to_dict(self):
        return {**self.stream.data,
                **self.unit.data,
                **self.system.data}
    
    @classmethod
    def from_flowsheets(cls, ID, flowsheets):
        """Return a new flowsheet with all registered objects from the given flowsheets."""
        new = cls(ID)
        isa = isinstance
        for flowsheet in flowsheets:
            if isa(flowsheet, str):
                flowsheet = cls.flowsheet[flowsheet]
            new.update(flowsheet)
        return new
    
    def diagram(self, kind: Optional[int|str]=None, file: Optional[str]=None, 
                format: Optional[str]=None, display: Optional[bool]=True,
                number: Optional[bool]=None, profile: Optional[bool]=None,
                label: Optional[bool]=None, title: Optional[str]=None,
                **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of
        all unit operations.

        Parameters
        ----------
        kind :
            * 0 or 'cluster': Display all units clustered by system.
            * 1 or 'thorough': Display every unit within the path.
            * 2 or 'surface': Display only elements listed in the path.
            * 3 or 'minimal': Display a single box representing all units.
        file : 
            File name to save diagram. 
        format:
            File format (e.g. "png", "svg"). Defaults to 'png'
        display : 
            Whether to display diagram in console or to return the graphviz
            object.
        number : 
            Whether to number unit operations according to their
            order in the system path.
        profile : 
            Whether to clock the simulation time of unit operations.
        label : 
            Whether to label the ID of streams with sources and sinks.

        """
        if title is None: title = ''
        return self.create_system(None).diagram(kind or 'thorough', file, format,
                                                display, number, profile, label,
                                                title, **graph_attrs)
    
    def create_system(self, ID: Optional[str]="", 
                      ends: Optional[Iterable[AbstractStream]]=None,
                      operating_hours: Optional[float]=None,
                      **kwargs):
        """
        Create a System object from all units and streams defined in the flowsheet.
        
        Parameters
        ----------
        ID : 
            Name of system.
        ends : 
            End streams of the system which are not products. Specify this
            argument if only a section of the complete system is wanted, or if
            recycle streams should be ignored.
        operating_hours : 
            Number of operating hours in a year. This parameter is used to
            compute annualized properties such as utility cost and material cost
            on a per year basis.
        
        """
        return System.from_units(ID, self.unit, ends, 
                                 operating_hours, **kwargs)
    
    def __call__(self, ID: str|type[AbstractUnit], strict: Optional[bool]=False):
        """
		Return requested biosteam item or a list of all matching items.
    
        Parameters
        ----------
        ID :
            ID of the requested item or Unit subclass.
        strict : 
            Whether an exact match is required. 
            
        """
        isa = isinstance
        if isa(ID, str):
            ID = ID.replace(' ', '_')
            obj = (self.stream.search(ID)
                   or self.unit.search(ID)
                   or self.system.search(ID))
            if not obj:
                if strict:
                    raise LookupError(f"no registered item '{ID}'")
                else:
                    obj = [i for i in self.unit if ID in ' '.join([i.__class__.__name__, i.ID])]
                    N = len(obj)
                    if N == 0:
                        raise LookupError(f"no registered item '{ID}'")
                    elif N == 1:
                        obj = obj[0]
            return obj
        elif issubclass(ID, bst.Unit):
            cls = ID
            obj = [i for i in self.unit if isa(i, cls)]
            N = len(obj)
            if N == 0:
                raise LookupError(f"no registered item '{ID}'")
            elif N == 1:
                obj = obj[0]
            else:
                obj = sorted(obj, key=lambda x: x.ID)
            return obj
        else:
            raise TypeError('ID must be either a string or a Unit subclass')
    
    def __str__(self):
        return self.ID
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'


class MainFlowsheet(Flowsheet):
    """
	Create a MainFlowsheet object which automatically registers 
    biosteam objects as they are created. For a tutorial on flowsheets,
    visit :doc:`tutorial/Managing_flowsheets`.
	"""
    
    __slots__ = ()
    
    line = "Main flowsheet"
        
    def set_flowsheet(self, flowsheet, new=False):
        """Set main flowsheet that is updated with new biosteam objects."""
        if isinstance(flowsheet, Flowsheet):
            dct = flowsheet.__dict__
        elif isinstance(flowsheet, str):
            if not new and flowsheet in self.flowsheet:
                dct = main_flowsheet.flowsheet[flowsheet].__dict__
            else:
                new_flowsheet = Flowsheet(flowsheet)
                self.flowsheet.__dict__[flowsheet] = new_flowsheet
                dct = new_flowsheet.__dict__
        else:
            raise TypeError('flowsheet must be a Flowsheet object')
        AbstractStream.registry = dct['stream']
        System.registry = dct['system']
        AbstractUnit.registry = dct['unit']
        object.__setattr__(self, '__dict__', dct)
        
    def get_flowsheet(self):
        return self.flowsheet[self.ID]
    
    def __new__(cls, ID):
        main_flowsheet.set_flowsheet(ID)
        return main_flowsheet

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'
    
    
#: Main flowsheet where objects are registered by ID.
#: Use the `set_flowsheet` to change the main flowsheet.
F = main_flowsheet = object.__new__(MainFlowsheet)
main_flowsheet.set_flowsheet(
    Flowsheet.from_registries(
        'default', AbstractStream.registry, AbstractUnit.registry, System.registry
    )
)
