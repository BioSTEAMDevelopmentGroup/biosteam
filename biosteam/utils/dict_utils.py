# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:03:08 2018

This module includes classes and functions relating dictionaries.

@author: Guest Group
"""
from biosteam.exceptions import IDconflict
from weakref import ref

__all__ = ('getDict', 'WeakRefBook', 'get_doc_units', 'merge', 'get_vals')


# %% Dictionary Classes

class getDict(dict):
    """Create a dictionary that return None instead of raising a KeyError."""

    def __getitem__(self, key):
        return self.get(key)

    def __repr__(self):
        return type(self).__name__ + '\n' + super().__repr__().replace('], ', '],\n ')


class WeakRefBook(dict):
    """Create a WeakRefBook object that stores only weak references. When item getter is used, it returns the call of the weak reference. If the key does not exist, it returns None. When check_ID is True, an IDconflict error will be raised when setting with a key/ID already occupied by an object.

    **Parameters**

         **check_ID:** [bool] if True, no key/ID can be specified twice.

         **kwargs: [dict] Key-value pairs to initialize dictionary.

    """

    __slots__ = ['name', 'check_ID']

    def __init__(self, check_ID=False, **kwargs):
        self.check_ID = check_ID
        for key, item in kwargs:
            self[key] = item

    def __setitem__(self, ID, obj):
        if self.check_ID and self[ID]:
            raise IDconflict(f"A '{type(obj).__name__}' object with key/ID, '{ID}', already exists. Add a '*' to ID to ignore possible conflicts.")
        if type(obj) is not ref:
            obj = ref(obj)
        super().__setitem__(ID, obj)

    def __delitem__(self, key):
        if key in self:
            super().__delitem__(key)

    def __getitem__(self, key):
        return self.get(key, lambda: None)()

    def __repr__(self):
        return f'{type(self).__name__}{tuple(self.keys())}'


# %% Dictionary functions

def get_doc_units(udict, doc):
    """Update dictionary of units of measure from docstring.
    
    **Parameters**
    
        **doc:** [str] Keys must be in single quotes, '', and units of measure must follow in parentesis, ().

    **Returns**
    
        udict: [dict] Units of meaure.

    Example:
        
        .. code-block:: python
        
           >>> get_units({},
           ...           '''
           ...           'Duty': (kJ/hr)
           ...           'T': Temperature (K) 
           ...           ''')
           {'Duty': 'kJ/hr', 'T': 'K'}, {'T': 'Temperature'}
        
    """
    if not doc:
        return udict
    doc_list = doc.split("'")
    del doc_list[0]
    for i in range(int(len(doc_list)/2)):
        i = i*2
        key = doc_list[i]
        val = doc_list[i+1]
        
        # Get units
        par1 = val.find('(')
        par2 = val.find(')')
        if par1==-1 or par2==-1:
            unit = ''
        else:
            unit = val[par1+1:par2]
        if unit:
            udict[key] = unit        
        

def merge(*dict_args: dict) -> 'Merged dictionary':
    """Given any number of dictionaries, shallow copy and merge into a new dictionary, precedence goes to key value pairs in latter dicts."""
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def get_vals(dct: dict, *keys) -> 'generator(dict values)':
    """Get all key values as a generator."""
    return (dct[key] for key in keys)





