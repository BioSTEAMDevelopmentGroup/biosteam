# -*- coding: utf-8 -*-
"""
As BioSTEAM objects are created, weakrefs are automatically added to dictionaries. When find is called, it looks up the item in the dictionaries and returns it. The find function allows the user to find any BioSTEAM object during debbuging or at the console level.
"""
from biosteam.utils import WeakRefBook, getDict


# %% Flowsheet search

def find(item_ID) -> 'item':
    """Return requested biosteam item.

    **Parameters**

         **item_ID:** [str] ID of the requested item.

    """
    for dct in dictionaries:
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

dictionaries = (find.system,
                find.line,
                find.unit,
                find.stream)

