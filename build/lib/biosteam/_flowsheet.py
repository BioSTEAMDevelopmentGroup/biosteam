# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. The `find` object allows the user to find any Unit, Stream or System instance.  When `find` is called, it simply looks up the item and returns it. 

"""
from graphviz import Digraph
from IPython import display
from ._utils import Register, search_register

__all__ = ('find', 'Flowsheet')


def make_digraph(units, streams):
    """Return digraph of units and streams."""
    # Create a digraph and set direction left to right
    f = Digraph(format='svg')
    f.attr(rankdir='LR')
    # Set up unit nodes
    UD = {}  # Contains full description (ID and line) by ID
    for u in units:
        if hasattr(u, 'link_streams'): u.link_streams()
        graphics = u._graphics
        if not graphics.in_system:
            continue  # Ignore Unit
        
        # Initialize graphics and make Unit node with attributes
        Type = graphics.node_function(u) or u.line
        name = u.ID + '\n' + Type
        f.attr('node', **u._graphics.node)
        f.node(name)
        UD[u] = name
        
    keys = UD.keys()

    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', splines='normal', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')
    
    for s in streams:
        if not s: continue  # Ignore stream

        oU = s._source
        if oU:
            oi = oU._outs.index(s) 
        
        dU = s._sink
        if dU:
            di = dU._ins.index(s)
        
        # Make stream nodes / unit-stream edges / unit-unit edges
        if oU not in keys and dU not in keys: pass
            # Stream is not attached to anything
        elif oU not in keys:
            # Feed stream case
            f.node(s.ID)
            edge_in = dU._graphics.edge_in
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(s.ID, UD[dU])
        elif dU not in keys:
            # Product stream case
            f.node(s.ID)
            edge_out = oU._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(UD[oU], s.ID)
        else:
            # Process stream case
            edge_in = dU._graphics.edge_in
            edge_out = oU._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='normal',
                   **edge_in[di], **edge_out[oi])
            f.edge(UD[oU], UD[dU], label=s.ID)
    return f

def save_digraph(digraph, file, format):
    if not file:
        if format == 'svg':
            x = display.SVG(digraph.pipe(format=format))
        else:
            x = display.Image(digraph.pipe(format='png'))
        display.display(x)
    else:
        if '.' not in file:
            file += '.' + format
        img = digraph.pipe(format=format)
        f = open(file, 'wb')
        f.write(img)
        f.close()


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
    """Create a Flowsheet object which stores references to all stream, unit, and system objects."""
    
    #: [Register] All flowsheets.
    flowsheet = Flowsheets()
    
    def __init__(self, ID):        
        #: [Register] Contains all System objects as attributes.
        self.system = Register()
        
        #: [Register] Contains all Unit objects as attributes.
        self.unit = Register()
        
        #: [Register] Contains all Stream objects as attributes.
        self.stream = Register()
        
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
    
    def diagram(self, kind='surface'):
        """Display all units and attached streams.
        
        **kind:** Must be one of the following:
                * **'thorough':** Thoroughly display every unit.
                * **'surface':** Display units and recycle systems.
                * **'minimal':** Minimally display flowsheet as a box with feeds and products.
        
        """
        if kind == 'thorough':
            return self._thorough_diagram()
        elif kind == 'surface':
            return self._surface_diagram()
        elif kind == 'minimal':
            return self._minimal_diagram()
        else:
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'.")
    
    def _thorough_diagram(self, file=None, format='svg'):
        units = list(self.unit)
        units.reverse()
        streams = set()
        for u in units:
            streams.update(u._ins)
            streams.update(u._outs)
        f = make_digraph(units, streams)
        save_digraph(f, file, format)
    
    def _minimal_diagram(self):
        from . import _system, Stream
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
        unit.diagram(1)
        
    def _surface_diagram(self, file=None, format='svg'):
        from . import _system, Stream
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
                        refresh_units.add(source)
                    elif sink in i.units and source not in i.units:
                        if source: ins.append(s)
                        else: feeds.append(s)
                        refresh_units.add(sink)
                
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
        sys._thorough_diagram(file, format)
        for i in refresh_units:
            i._ins[:] = i._ins
            i._outs[:] = i._outs
    
    def __call__(self, ID):
        """Return requested biosteam item.
    
        **Parameters**
    
             **ID:** [str] ID of the requested item.
    
        """
        ID = ID.replace(' ', '_')
        obj = (search_register(self.stream, ID)
               or search_register(self.unit, ID)
               or search_register(self.system, ID))
        if not obj: raise LookupError(f"no registered item '{ID}'")
        return obj
    
    def __str__(self):
        return self.ID
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'


class MainFlowsheet(Flowsheet):
    """Create a MainFlowsheet object which automatically registers biosteam objects as they are created."""
    __slots__ = ()
    
    @staticmethod
    def set_flowsheet(flowsheet):
        """Set main flowsheet that is updated with new biosteam objects."""
        if isinstance(flowsheet, Flowsheet):
            object.__setattr__(find, '__dict__', flowsheet.__dict__)
        else:
            raise TypeError('flowsheet must be a Flowsheet object')
    
    __setattr__ = Flowsheets.__setattr__
    
    def __new__(cls):
        raise TypeError('cannot create new find object.')

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'
    
    
#: [find] Find BioSTEAM objects by ID.
find = object.__new__(MainFlowsheet)
find.set_flowsheet(Flowsheet('Default'))
    

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