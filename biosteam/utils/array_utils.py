# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:42:49 2018

This module includes classes and functions relating arrays.

@author: Guest Group
"""
import numpy as np
from array_collections import array, tuple_array, material_array, \
                              property_array, same_type_tuple, \
                              organized_list, same_type, reorder, get_frac
from free_properties import PropertyFactory

__all__ = ('array', 'tuple_array', 'material_array', 'property_array',
           'PropertyFactory', 'same_type_tuple', 'organized_list',
           'same_type', 'reorder', 'get_frac')

ndarray = np.ndarray
asarray = np.asarray

material_array._units = 'kmol/hr'