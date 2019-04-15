# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. The `find` object allows the user to find any Unit, Stream or System instance.  When `find` is called, it simply looks up the item and returns it. 

:doc:`Find unit operations and manage flowsheets` 

"""
from graphviz import Digraph
# import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
from IPython import display
# import threading
# import time
# import os

__all__ = ('find', 'stream_connector', 'Flowsheet')


def make_digraph(units, streams=None):
    """Return digraph of units and streams."""
    # Create a digraph and set direction left to right
    f = Digraph(format='svg')
    f.attr(rankdir='LR')
    # Set up unit nodes
    U = {}  # Contains units by ID
    UD = {}  # Contains full description (ID and line) by ID
    for u in units:
        u._link_streams()
        graphics = u._graphics
        if not graphics.in_system:
            continue  # Ignore Unit
        
        # Initialize graphics and make Unit node with attributes
        Type = graphics.node_function(u) or u.line
        name = u.ID + '\n' + Type
        f.attr('node', **u._graphics.node)
        f.node(name)
        U[u.ID] = u
        UD[u.ID] = name
        
    keys = UD.keys()

    # Set attributes for graph and streams
    f.attr('node', shape='rarrow', fillcolor='#79dae8',
           style='filled', orientation='0', width='0.6',
           height='0.6', color='black', peripheries='1')
    f.attr('graph', splines='normal', overlap='orthoyx',
           outputorder='edgesfirst', nodesep='0.15', maxiter='1000000')
    f.attr('edge', dir='foward')

    connections = set()
    if not streams:
        streams = set()
        for u in units:
            streams.update(u._ins)
            streams.update(u._outs)
        streams.difference_update(find._upstream_connections)
    for s in streams:
        if s.ID == 'Missing Stream':
            continue  # Ignore stream
        
        oU = s._source
        if oU:
            try: oi = oU._outs.index(s) 
            except:
                unit_connection = (oU, s)
                if unit_connection not in connections:
                    oi = oU._outs.index(s._upstream_connection) # The stream says it's source is that unit, but it really means that the unit is the source of the connection stream
                    connections.add(unit_connection)
            oU = oU._ID
        else: oi = None
        
        dU = s._sink
        if dU:
            try:
                di = dU._ins.index(s)
                dU = dU._ID
            except: # The stream is the upstream connection stream
                unit_connection = (dU, s)
                if unit_connection not in connections:
                    di = dU._ins.index(s._downstream_connection) # The stream says it's source is that unit, but it really means that the unit is the source of the connection stream
                    connections.add(unit_connection)
        else: di = None

        # Make stream nodes / unit-stream edges / unit-unit edges
        if oU not in keys and dU not in keys: pass
            # Stream is not attached to anything
        elif oU not in keys:
            # Feed stream case
            f.node(s.ID)
            edge_in = U[dU]._graphics.edge_in
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(s.ID, UD[dU])
        elif dU not in keys:
            # Product stream case
            f.node(s.ID)
            edge_out = U[oU]._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(UD[oU], s.ID)
        else:
            # Process stream case
            edge_in = U[dU]._graphics.edge_in
            edge_out = U[oU]._graphics.edge_out
            f.attr('edge', arrowtail='none', arrowhead='normal',
                   **edge_in[di], **edge_out[oi])
            f.edge(UD[oU], UD[dU], label=s.ID)
    return f


# %% Flowsheet search      
    
class Flowsheet:
    """Create a Flowsheet object which stores references to all stream, unit, and system objects."""
    #: [dict] All flowsheets
    flowsheet = {}
    
    def __init__(self, ID):
        #: [str] ID of flowsheet
        self.ID = ID
        
        #: [dict] Dictionary of systems
        self.system = {}
        
        #: [dict] Dictionary of units
        self.unit = {}
        
        #: [dict] Dictionary of streams
        self.stream = {}
        
        self.flowsheet[ID] = self
    
    def diagram(self):
        """Display all units and attached streams."""
        units = list(self.unit.values())
        units.reverse()
        img = make_digraph(units).pipe('png')
        display.display(display.Image(img))
    
    def __call__(self, item_ID) -> 'item':
        """Return requested biosteam item.
    
        **Parameters**
    
             **item_ID:** [str] ID of the requested item.
    
        """
        item_ID = item_ID.replace(' ', '_')
        obj = (self.stream.get(item_ID)
               or self.unit.get(item_ID)
               or self.system.get(item_ID))
        if not obj: raise ValueError(f"No registered item '{item_ID}'")
        return obj
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'


class find(Flowsheet):
    """Create a find object which can search through flowsheets."""
    __slots__ = ()
    _upstream_connections = set()
    
    @property
    def mainflowsheet(self):
        """[Flowsheet] Main flowsheet that is updated with new biosteam objects"""
        return find._mainflowsheet
    
    @mainflowsheet.setter
    def mainflowsheet(self, flowsheet):
        if isinstance(flowsheet, Flowsheet):
            find.__dict__ = flowsheet.__dict__
        else:
            raise TypeError('Main flowsheet must be a Flowsheet object')
        find._mainflowsheet = flowsheet
    
    def __new__(cls):
        raise TypeError('Cannot create new Find object.')

    def __repr__(self):
        return f'<{type(self).__name__}: mainflowsheet={self.ID}>'
    
    
#: [find] Find BioSTEAM objects by ID.
find = object.__new__(find)
find.mainflowsheet = Flowsheet('Default')


# %% Connect between different flowsheets

def stream_connector(upstream, downstream):
    """Return a function that copies specifications from `upstream` to `downstream`. This serves to connect different flowsheets.
    
    **Parameters**
    
        **upstream:** [Stream] Non-zero species of this stream should be specified in the species object of `downstream`.
        
        **downstream:** [Stream] Flow rate, T, P, and phase information will be copied from `upstream` to this stream.
    
    """
    # Source and sink.
    upstream._sink = downstream._sink
    downstream._source = upstream._source
    downstream._upstream_connection = upstream
    upstream._downstream_connection = downstream
    find._upstream_connections.add(upstream)
    def connect():
        # Flow rate, T, P and phase
        index, species = upstream.nonzero_species
        downstream.setflow(upstream._molarray[index], species)
        downstream.T = upstream.T
        downstream.P = upstream.P
        downstream.phase = upstream.phase
    return connect


        
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