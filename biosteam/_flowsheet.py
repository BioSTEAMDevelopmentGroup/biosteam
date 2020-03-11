# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. The `find` object allows the user to find any Unit, Stream or System instance.  When `find` is called, it simply looks up the item and returns it. 

"""
from thermosteam.utils import Registry
from ._digraph import make_digraph, save_digraph
from thermosteam import Stream
from ._unit import Unit
from ._facility import Facility
from ._system import System

__all__ = ('find', 'Flowsheet')

# %% Flowsheet search      

class Flowsheets:
    __getattribute__ = __getitem__ = object.__getattribute__
    
    def __setattr__(self, key, value):
        raise TypeError(f"'{type(self).__name__}' object does not support attribute assignment")
    __setitem__ = __setattr__
    
    def __iter__(self):
        yield from self.__dict__.values()
    
    def __delattr__(self, key):
        if key == find.ID:
            raise AttributeError('cannot delete main flowsheet')
        else:
            super().__delattr__(key)
    
    def __repr__(self):
        return f'Register:\n ' + '\n '.join([repr(i) for i in self])
    
    
class Flowsheet:
    """
    Create a Flowsheet object which stores references to all stream, unit,
    and system objects."""
    
    #: [Register] All flowsheets.
    flowsheet = Flowsheets()
    
    def __init__(self, ID):        
        #: [Register] Contains all System objects as attributes.
        self.system = Registry()
        
        #: [Register] Contains all Unit objects as attributes.
        self.unit = Registry()
        
        #: [Register] Contains all Stream objects as attributes.
        self.stream = Registry()
        
        #: [str] ID of flowsheet.
        self._ID = ID
        self.flowsheet.__dict__[ID] = self
    
    def __setattr__(self, key, value):
        if hasattr(self, '_ID'):
            raise TypeError(f"'{type(self).__name__}' object does not support attribute assignment")
        else:
            super().__setattr__(key, value)
    
    @property
    def ID(self):
        return self._ID
    
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
        _system._streamUnit('\n'.join([i.ID for i in feeds]),
                           None, feed)
        _system._streamUnit('\n'.join([i.ID for i in products]),
                           product, None)
        unit = _system._systemUnit(self.ID, feed, product)
        unit.line = 'flowsheet'
        unit.diagram(1, file, format, **graph_attrs)
        
    def _surface_diagram(self, file, format, **graph_attrs):
        from . import _system
        units = set(self.unit)
        StrUnit = _system._streamUnit
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
                
                subsystem_unit = _system._systemUnit(i.ID, ins, outs)
                units.difference_update(i.units)
                units.add(subsystem_unit)
        
        sys = _system.System(None, units)
        sys._thorough_diagram(file, format, **graph_attrs)
        for u, ins, outs in refresh_units:
            u._ins[:] = ins
            u._outs[:] = outs
    
    def create_system(self, ID=None, facilities=(), end_streams=()):
        streams = tuple(self.stream)
        feeds = [i for i in streams if i._sink and not 
                 (i._source or isinstance(i._sink, Facility))]
        feeds_negflows = [(i, -i.F_mass) for i in feeds]
        feeds_negflows = sorted(feeds_negflows, key=lambda x: x[1])
        feedstock, *feeds = [i[0] for i in feeds_negflows]
        return System.from_feedstock(ID or self.ID + '_sys', feedstock, feeds,
                                     facilities, end_streams)
    
    def __call__(self, ID):
        """Return requested biosteam item.
    
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
    """Create a MainFlowsheet object which automatically registers 
    biosteam objects as they are created."""
    __slots__ = ()
    
    @staticmethod
    def set_flowsheet(flowsheet):
        """Set main flowsheet that is updated with new biosteam objects."""
        if isinstance(flowsheet, Flowsheet):
            dct = flowsheet.__dict__
        elif isinstance(flowsheet, str):
            if flowsheet in find.flowsheet:
                dct = find.flowsheet[flowsheet].__dict__
            else:
                new_flowsheet = Flowsheet(flowsheet)
                find.flowsheet.__dict__[flowsheet] = new_flowsheet
                dct = new_flowsheet.__dict__
        else:
            raise TypeError('flowsheet must be a Flowsheet object')
        Stream.registry = dct['stream']
        System.registry = dct['system']
        Unit.registry = dct['unit']
        object.__setattr__(find, '__dict__', dct)
        
    __setattr__ = Flowsheets.__setattr__
    
    def __new__(cls):
        raise TypeError('cannot create new find object.')

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'
    
    
#: [find] Find BioSTEAM objects by ID.
find = object.__new__(MainFlowsheet)
find.set_flowsheet('default')
    

# %% Attempt at contantly rendered digraph

# #: [bool] True if flowsheet is contantly rendered
# self.live = True

# def _show_once(self):
    #     fig = plt.figure()
    #     plt.axis('off')
    #     f = open('diagram.png', 'wb')
    #     diagram = make_digraph(self.unit.values()).pipe(format='png')
    #     f.write(diagram)
    #     img = mpimg.imread('diagram.png')
    #     plt.imshow(img, interpolation='spline36')
    #     fig.show()
    #     # img = Image.open('diagram.png')
    #     # img.show() 
    
    # def _show(self):
    #     fig = plt.figure()
    #     plt.axis('off')
    #     f = open('diagram.png', 'wb')
    #     for i in range(3):
    #         diagram = make_digraph(self.unit.values()).pipe(format='png')
    #         f.write(diagram)
    #         img = mpimg.imread('diagram.png')
    #         plt.clf() 
    #         plt.imshow(img, interpolation='spline36')
    #         fig.show()
    #         time.sleep(10)
    #     f.close()
    #     os.remove('diagram.png')
    
    # def show(self):
    #     t1 = threading.Thread(target=self._show())
    #     t1.start()