# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. The `main_flowsheet` object allows the user to find any Unit, Stream or System instance.  When `main_flowsheet` is called, it simply looks up the item and returns it. 
"""
import sys
from PyQt5.QtWidgets import QApplication
from thermosteam.utils import Registry
from ._viz import FlowsheetWidget, make_digraph, save_digraph
from thermosteam import Stream
from ._unit import Unit
from ._facility import Facility
from ._system import System
from ._network import Network
from PyQt5 import QtCore

__all__ = ('main_flowsheet', 'Flowsheet')

# %% Functions

def get_feeds_big_to_small(streams):
    isa = isinstance
    feeds = [i for i in streams if i._sink and not 
             (i._source or isa(i._sink, Facility))]
    feeds_negflows = [(i, -i.F_mass) for i in feeds]
    feeds_negflows = sorted(feeds_negflows, key=lambda x: x[1])
    return [i[0] for i in feeds_negflows]

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
    :doc:`tutorial/Managing flowsheets.ipynb`.
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
            return self._thorough_diagram(file, format, **graph_attrs)
        elif kind == 'surface':
            return self._surface_diagram(file, format, **graph_attrs)
        elif kind == 'minimal':
            return self._minimal_diagram(file, format, **graph_attrs)
        else:
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'.")
    
    def _thorough_diagram(self, file, format, **graph_attrs):
        units = list(self.unit)
        units.reverse()
        streams = set()
        for u in units:
            streams.update(u._ins)
            streams.update(u._outs)
        f = make_digraph(units, streams, format=format, **graph_attrs)
        save_digraph(f, file, format)
    
    def _minimal_diagram(self, file, format, **graph_attrs):
        from . import _system
        streams = list(self.stream)
        feeds = set(filter(_system._isfeed, streams))
        products = set(filter(_system._isproduct, streams))
        product = Stream(None)
        product._ID = ''
        feed = Stream(None)
        feed._ID = ''
        _system.StreamUnit('\n'.join([i.ID for i in feeds]),
                           None, feed)
        _system.StreamUnit('\n'.join([i.ID for i in products]),
                           product, None)
        unit = _system.SystemUnit(self.ID, feed, product)
        unit.line = 'flowsheet'
        unit.diagram(1, file, format, **graph_attrs)
        
    def _surface_diagram(self, file, format, **graph_attrs):
        from . import _system
        units = set(self.unit)
        StrUnit = _system.StreamUnit
        refresh_units = set()
        for i in self.system:
            if i.recycle and not any(sub.recycle for sub in i.subsystems):
                outs = []
                ins = []
                feeds = []
                products = []
                for s in i.streams:
                    source = s._source
                    sink = s._sink
                    if source in i.units and sink not in i.units:
                        if sink: outs.append(s)
                        else: products.append(s)
                        u_io = (source, tuple(source.ins), tuple(source.outs))
                        refresh_units.add(u_io)
                    elif sink in i.units and source not in i.units:
                        if source: ins.append(s)
                        else: feeds.append(s)
                        u_io = (sink, tuple(sink.ins), tuple(sink.outs))
                        refresh_units.add(u_io)
                
                if len(feeds) > 1:
                    feed = Stream(None)
                    feed._ID = ''
                    units.add(StrUnit('\n'.join([i.ID for i in feeds]),
                                      None, feed))
                    ins.append(feed)
                else: ins += feeds
                
                if len(products) > 1:
                    product = Stream(None)
                    product._ID = ''
                    units.add(StrUnit('\n'.join([i.ID for i in products]),
                                      product, None))
                    outs.append(product)
                else: outs += products
                
                subsystem_unit = _system.SystemUnit(i.ID, ins, outs)
                units.difference_update(i.units)
                units.add(subsystem_unit)
        
        sys = _system.System(None, units)
        sys._thorough_diagram(file, format, **graph_attrs)
        for u, ins, outs in refresh_units:
            u._ins[:] = ins
            u._outs[:] = outs
    
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
        feedstock, *feeds = get_feeds_big_to_small(feeds or self.stream)
        facilities = self.get_facilities()
        return System.from_feedstock(ID, feedstock, feeds,
                                     facilities, ends)
    
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
        feedstock, *feeds = get_feeds_big_to_small(feeds or self.stream)
        return Network.from_feedstock(feedstock, feeds, ends)
    
    def create_path(self, feeds=None, ends=()):
        isa = isinstance
        network = self.create_network(feeds, ends)
        net2sys = System.from_network
        return tuple([(net2sys('', i) if isa(i, Network) else i)
                      for i in network.path])
    
    def get_facilities(self):
        isa = isinstance
        return [i for i in self.unit if isa(i, Facility)]
    
    def view(self, autorefresh=True):
        """Create an interactive process flowsheet diagram that autorefreshes itself."""
        widget = FlowsheetWidget(self, autorefresh)
        widget.show()
        return widget
    
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
    visit :doc:`tutorial/Managing flowsheets.ipynb`.
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
