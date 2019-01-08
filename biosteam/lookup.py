# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:20:32 2018

@author: yoelr
"""
from biosteam.utils import WeakRefBook, getDict


# %% Flowsheet search
class metaLookUp(type):
    """Meta class for defining the LookUp call."""
    def __call__(cls, ID: 'item ID') -> 'item':
        """Return requested item from LookUp's weakref dictionaries.

        **Parameters**

             **ID:** [str] ID of the requested item.

        """
        for dct_name in ('system', 'line', 'unit', 'stream'):
            obj = getattr(cls, dct_name)[ID]
            if obj:
                return obj

    @property
    def check_ID(cls):
        return cls._check_ID

    @check_ID.setter
    def check_ID(cls, boolean):
        cls._check_ID = boolean
        for dct_name in ('system', 'line', 'unit', 'stream'):
            getattr(cls, dct_name).check_ID = boolean


class LookUp(metaclass=metaLookUp):
    """This class stores weekrefs for all systems, units, streams and lines (Unit classes) in class attribute dictionaries. As BioSTEAM objects are created, weakrefs are automatically added to the dictionaries. When called, it looks up the item in the dictionaries and returns it.

    **Parameters**

        **ID:** [str] ID of the requested item.

    **Return**

        Requested item or None.

    **Class Attributes**

        check_ID = False: [bool] if True, no object ID can be specified twice, preventing objects with the same name. If false, two or more objects can exist with the same ID.

        stream = {}: [WeakRefBook] Dictionary of Stream objects.

        unit = {}: [WeakRefBook] Dictionary of Unit objects.

        system = {}: [WeakRefBook] Dictionary of System objects.

        line = {}: [getDict] Dictionary of lines and Unit classes (each line contains a list of Unit classes that inhearit from the same line).

    .. Note::
        Does not create a new object when called. Simply returns requested item.

     """
    # [bool] if True, no object ID can be specified twice, preventing objects with the same name. If false, two or more objects can exist with the same ID.
    _check_ID = False  
    
    # [WeakRefBook] Dictionary of systems
    system = WeakRefBook()  
    
    # [getDict] Dictionary of lines (each line contains a list of Unit classes that inhearit from the same line)
    line = getDict()        
    
    # [WeakRefBook] Dictionary of units
    unit = WeakRefBook()    
    
    # [WeakRefBook] Dictionary of streams
    stream = WeakRefBook()  

    # Filler init function for Sphinx documentation
    def __init__(cls, ID: 'item ID') -> 'item':
        """Return requested item from LookUp's weakref dictionaries.

        **Parameters**

             ID: [str] ID of the requested item.

        """

    @classmethod
    def empty(cls, dct_name='all'):
        """Empty the specified dictionary. This prevents ID conflicts.

        **Parameters**

             dct_name: [str] Name of dictionary to be emptied. Can be 'system', 'unit', 'stream', or 'all' for all dictionaries mentioned.

        """
        if dct_name == 'all':
            for dct_name in ('system', 'unit', 'stream'):
                cls.empty(dct_name)
        else:
            setattr(cls, dct_name, cls._WeakRefBook(dct_name))
        if dct_name == 'stream':
            from .Stream import Stream
            Stream._default_ID[1] = 0
