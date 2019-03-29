# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, weakrefs are automatically added to dictionaries. When find is called, it looks up the item in the dictionaries and returns it. The find function allows the user to find any BioSTEAM object during debbuging or at the console level.
"""
__all__ = ('find', 'stream_connector', 'Flowsheet')

# %% Flowsheet search

class metaFlowsheet(type):
    @property
    def main(self):
        """[Flowsheet] Main flowsheet that is updated with new biosteam objects"""
        return self._main
    
    @main.setter
    def main(self, flowsheet):
        if isinstance(flowsheet, Flowsheet):
            Flowsheet._main = flowsheet
        else:
            raise TypeError('main must be a Flowsheet object')
            
    
class Flowsheet(metaclass=metaFlowsheet):
    """Create a Flowsheet object which stores references to all stream, unit, and system objects."""
    # [Dict] All flowsheets
    flowsheets = {}
    def __init__(self, ID):
        #: [str] ID of flowshett
        self.ID = ID
        
        # [dict] Dictionary of systems
        self.system = sys = {}
        
        # [dict] Dictionary of lines (each line contains a list of Unit classes that inhearit from the same line)
        self.line = l = {}
        
        # [dict] Dictionary of units
        self.unit = u = {}
        
        # [dict] Dictionary of streams
        self.stream = s = {}
        
        #: All search dictionaries
        self._dicts = (s, u, sys, l)
        
        self.main = self
        
        self.flowsheets[ID] = self
    
    def __call__(self, item_ID) -> 'item':
        """Return requested biosteam item.
    
        **Parameters**
    
             **item_ID:** [str] ID of the requested item.
    
        """
        for dct in self._dicts:
            obj = dct.get(item_ID)
            if obj: return obj
        print(f"No registered item '{item_ID}'")
    
    main = metaFlowsheet.main
    
    @classmethod    
    def find(cls, item_ID) -> 'item':
        """Return requested biosteam item from main flowsheet.
    
        **Parameters**
    
             **item_ID:** [str] ID of the requested item.
    
        """
        return cls._main(item_ID.replace(' ', '_'))
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self.ID}>'

#: [Flowsheet] Default flowsheet
flowsheet = Flowsheet('default')
find = Flowsheet.find


# %% Connect between different flowsheets

def stream_connector(upstream, downstream):
    """Return a function that copies specifications from `upstream` to `downstream`. This serves to connect different flowsheets.
    
    **Parameters**
    
        **upstream:** [Stream] Non-zero species of this stream should be specified in the species object of `downstream`.
        
        **downstream:** [Stream] Flow rate, T, P, and phase information will be copied from `upstream` to this stream.
    
    """
    # Source and sink. Connection precedense goes to downstream
    upstream._sink = downstream._sink
    downstream._source = upstream._source
    downstream._upstream_connection = upstream
    def connect():
        # Flow rate, T, P and phase
        index, species = upstream.nonzero_species
        downstream.setflow(upstream._molarray[index], species)
        downstream.T = upstream.T
        downstream.P = upstream.P
        downstream.phase = upstream.phase
    return connect


