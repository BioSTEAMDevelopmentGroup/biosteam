# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, they are automatically registered. `find` is a MainFlowsheet object that allows the user to find any Unit, Stream or System instance.  When `find` is called, it simply looks up the item and returns it. 

**Examples**

    .. code-block:: python

       >>> from biosteam import find
       >>> import lipidcane # You must pip install lipidcane first
       >>> find
       <MainFlowsheet: Lipidcane>
       >>> find('Lipid_cane')
       <Stream: Lipid_cane>

"""
__all__ = ('find', 'stream_connector', 'MainFlowsheet', 'Flowsheet')

# %% Flowsheet search      
    
class Flowsheet:
    """Create a Flowsheet object which stores references to all stream, unit, and system objects."""
    
    def __init__(self, ID):
        #: [str] ID of flowsheet
        self.ID = ID
        
        #: [dict] Dictionary of systems
        self.system = {}
        
        #: [dict] Dictionary of units
        self.unit = {}
        
        #: [dict] Dictionary of streams
        self.stream = {}
    
    def __call__(self, item_ID) -> 'item':
        """Return requested biosteam item.
    
        **Parameters**
    
             **item_ID:** [str] ID of the requested item.
    
        """
        item_ID = item_ID.replace(' ', '_')
        obj = self.stream.get(item_ID)
        if obj: return obj
        obj = self.unit.get(item_ID)
        if obj: return obj
        obj = self.system.get(item_ID)
        if obj: return obj
        print(f"No registered item '{item_ID}'")
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'


class MainFlowsheet(Flowsheet):
    """Create a Flowsheet object which stores references to all stream, unit, and system objects."""
    __slots__ = ()
    
    @property
    def flowsheet(self):
        """[Flowsheet] Main flowsheet that is updated with new biosteam objects"""
        return MainFlowsheet._flowsheet
    
    @flowsheet.setter
    def flowsheet(self, flowsheet):
        if isinstance(flowsheet, Flowsheet):
            find.__dict__ = flowsheet.__dict__
        else:
            raise TypeError('Main flowsheet must be a Flowsheet object')
        MainFlowsheet._flowsheet = flowsheet
    
    def __new__(cls):
        raise TypeError('Cannot create new MainFlowsheet object. Only one main flowsheet can exist.')

    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'
    
    
#: [MainFlowsheet] Main flowsheet where Stream, Unit, and System objects are registered.
find = object.__new__(MainFlowsheet)
find.flowsheet = Flowsheet('Default')


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
    def connect():
        # Flow rate, T, P and phase
        index, species = upstream.nonzero_species
        downstream.setflow(upstream._molarray[index], species)
        downstream.T = upstream.T
        downstream.P = upstream.P
        downstream.phase = upstream.phase
    return connect


