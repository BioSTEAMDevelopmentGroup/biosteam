# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, weakrefs are automatically added to dictionaries. When find is called, it looks up the item in the dictionaries and returns it. The find function allows the user to find any BioSTEAM object during debbuging or at the console level.
"""
from .utils import WeakRefBook, getDict

__all__ = ('find', 'stream_connector')

# %% Flowsheet search

def find(item_ID) -> 'item':
    """Return requested biosteam item.

    **Parameters**

         **item_ID:** [str] ID of the requested item.

    """
    for dct in _dicts:
        obj = dct[item_ID]
        if obj: return obj
    print(f"No registered item '{item_ID}'")

# [WeakRefBook] Dictionary of systems
find.system = WeakRefBook()  

# [getDict] Dictionary of lines (each line contains a list of Unit classes that inhearit from the same line)
find.line = getDict()        

# [WeakRefBook] Dictionary of units
find.unit = WeakRefBook()    

# [WeakRefBook] Dictionary of streams
find.stream = WeakRefBook()  

_dicts = (find.system,
          find.line,
          find.unit,
          find.stream)

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


