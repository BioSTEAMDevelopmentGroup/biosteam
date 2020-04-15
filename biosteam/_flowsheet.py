# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. The `main_flowsheet` object allows the user to find any Unit, Stream or System instance.  When `main_flowsheet` is called, it simply looks up the item and returns it. 
"""
from thermosteam.utils import Registry
from ._digraph import (digraph_from_units,
                       digraph_from_units_and_streams, 
                       finalize_digraph,
                       minimal_digraph,
                       update_surface_units)
from thermosteam import Stream
from ._unit import Unit
from ._facility import Facility
from ._system import System
from ._network import Network

try:
	from ._digraph.widget import FlowsheetWidget
except:
	pass

__all__ = ('main_flowsheet', 'Flowsheet')

# %% Functions

def get_feeds_from_streams(streams):
    isa = isinstance
    return [i for i in streams if i._sink and not 
            (i._source or isa(i._sink, Facility))]

def sort_feeds_big_to_small(feeds):
    feeds.sort(key=lambda feed: -feed.F_mass)
    

# %% Flowsheet search      

class Flowsheets:
    __getattribute__ = __getitem__ = object.__getattribute__
    
    def __setattr__(self, key, value):
        raise TypeError(f"'{type(self).__name__}' object does not support attribute assignment")
    __setitem__ = __setattr__
    
    def __iter__(self):
        yield from self.__dict__.values()
    
    def __delattr__(self, key):
        if key == main_flowsheet.ID:
            raise AttributeError('cannot delete main flowsheet')
        else:
            super().__delattr__(key)
    
    def __repr__(self):
        return f'Register:\n ' + '\n '.join([repr(i) for i in self])
    
    
class Flowsheet:
    """
    Create a Flowsheet object which stores references to all stream, unit,
    and system objects. For a tutorial on flowsheets, visit
    :doc:`tutorial/Managing_flowsheets`.
	"""
    line = "Flowsheet"
    
    #: [Register] All flowsheets.
    flowsheet = Flowsheets()
    
    def __new__(cls, ID):        
        self = super().__new__(cls)
        
        #: [Register] Contains all System objects as attributes.
        self.system = Registry()
        
        #: [Register] Contains all Unit objects as attributes.
        self.unit = Registry()
        
        #: [Register] Contains all Stream objects as attributes.
        self.stream = Registry()
        
        #: [str] ID of flowsheet.
        self._ID = ID
        self.flowsheet.__dict__[ID] = self
        return self
    
    def __reduce__(self):
        return self.from_registries, self.registries
    
    def __setattr__(self, key, value):
        if hasattr(self, '_ID'):
            raise TypeError(f"'{type(self).__name__}' object does not support attribute assignment")
        else:
            super().__setattr__(key, value)
    
    def view(self):
        """
        Create an interactive process flowsheet diagram
        that autorefreshes itself.
        """
        widget = FlowsheetWidget(self)
        widget.show()
        return widget
    
    @property
    def ID(self):
        return self._ID
    
    @classmethod
    def from_registries(cls, stream, unit, system):
        flowsheet = super().__new__(cls)
        flowsheet.stream = stream
        flowsheet.unit = unit
        flowsheet.system = system
        return flowsheet
    
    @property
    def registries(self):
        return (self.stream, self.unit, self.system)
    
    def discard(self, ID):
        for registry in self.registries:
            if ID in registry: 
                registry.discard(ID)
                break
    
    def update(self, flowsheet):
        for registry, other_registry in zip(self.registries, flowsheet.registries):
            registry.__dict__.update(other_registry.__dict__)
    
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
    
    def diagram(self, kind='surface', file=None, format='svg', **graph_attrs):
        """Display all units and attached streams.
        
        Parameters
        ----------
        kind='surface' : Must be one of the following
            * **'thorough':** Thoroughly display every unit.
            * **'surface':** Display units and recycle systems.
            * **'minimal':** Minimally display flowsheet as a box with feeds and products.
        
        """
        if kind == 'thorough':
            f = digraph_from_units_and_streams(self.unit, self.stream, 
                                               format=format, **graph_attrs)
        elif kind == 'surface':
            f = self._surface_digraph(format, **graph_attrs)
        elif kind == 'minimal':
            f = minimal_digraph(self.ID, self.units, self.streams, **graph_attrs)
        else:
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'.")
        finalize_digraph(f, file, format)
    
    def _surface_digraph(self, format, **graph_attrs):
        surface_units = set(self.unit)
        old_unit_connections = set()
        for i in self.system:
            if i.recycle and not any(sub.recycle for sub in i.subsystems):
                surface_units.difference_update(i.units)
                update_surface_units(i.ID, i.streams, i.units,
                                     surface_units, old_unit_connections)
        
        f = digraph_from_units(surface_units)
        for u, ins, outs in old_unit_connections:
            u._ins[:] = ins
            u._outs[:] = outs
        return f
    
    def create_system(self, ID="", feeds=None, ends=()):
        """
        Create a System object from all units and streams defined in the flowsheet.
        
        Parameters
        ----------
        ID : str, optional
            Name of system.
        ends : Iterable[:class:`~thermosteam.Stream`]
            End streams of the system which are not products. Specify this argument
			if only a section of the system is wanted, or if recycle streams should be 
			ignored.
        
        """
        feeds = get_feeds_from_streams(self.stream)
        if feeds:
            sort_feeds_big_to_small(feeds)
            feedstock, *feeds = feeds
            facilities = self.get_facilities()
            system = System.from_feedstock(ID, feedstock, feeds,
                                           facilities, ends)
        else:
            system = System(ID, ())
        return system
    
    def create_network(self, feeds=None, ends=()):
        """
        Create a Network object from all units and streams defined in the flowsheet.
        
        Parameters
        ----------
        ends : Iterable[:class:`~thermosteam.Stream`]
            End streams of the system which are not products. Specify this argument
			if only a section of the system is wanted, or if recycle streams should be 
			ignored.
        
        """
        feeds = get_feeds_from_streams(self.stream)
        if feeds:
            sort_feeds_big_to_small(feeds)
            feedstock, *feeds = feeds
            network = Network.from_feedstock(feedstock, feeds, ends)
        else:
            network = Network([])
        return network
    
    def create_path(self, feeds=None, ends=()):
        isa = isinstance
        network = self.create_network(feeds, ends)
        net2sys = System.from_network
        return tuple([(net2sys('', i) if isa(i, Network) else i)
                      for i in network.path])
    
    def get_facilities(self):
        isa = isinstance
        return [i for i in self.unit if isa(i, Facility)]
    
    def __call__(self, ID):
        """
		Return requested biosteam item.
    
        Parameters
        ----------
        ID : str
              ID of the requested item.
    
        """
        ID = ID.replace(' ', '_')
        obj = (self.stream.search(ID)
               or self.unit.search(ID)
               or self.system.search(ID))
        if not obj: raise LookupError(f"no registered item '{ID}'")
        return obj
    
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
    line = "Main flowsheet"
        
    def set_flowsheet(self, flowsheet):
        """Set main flowsheet that is updated with new biosteam objects."""
        if isinstance(flowsheet, Flowsheet):
            dct = flowsheet.__dict__
        elif isinstance(flowsheet, str):
            if flowsheet in self.flowsheet:
                dct = main_flowsheet.flowsheet[flowsheet].__dict__
            else:
                new_flowsheet = Flowsheet(flowsheet)
                self.flowsheet.__dict__[flowsheet] = new_flowsheet
                dct = new_flowsheet.__dict__
        else:
            raise TypeError('flowsheet must be a Flowsheet object')
        Stream.registry = dct['stream']
        System.registry = dct['system']
        Unit.registry = dct['unit']
        object.__setattr__(self, '__dict__', dct)
        
    def get_flowsheet(self):
        return self.flowsheet[self.ID]
        
    __setattr__ = Flowsheets.__setattr__
    
    def __new__(cls, ID):
        main_flowsheet.set_flowsheet(ID)
        return main_flowsheet

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'
    
    
#: [main_flowsheet] Main flowsheet where objects are registered by ID.
main_flowsheet = object.__new__(MainFlowsheet)
main_flowsheet.set_flowsheet('default')
