# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:03:08 2018

This module includes classes and functions relating dictionaries.

@author: Guest Group
"""
from biosteam.exceptions import IDconflict
from weakref import ref

__all__ = ('get_doc_units',)


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
        
        
