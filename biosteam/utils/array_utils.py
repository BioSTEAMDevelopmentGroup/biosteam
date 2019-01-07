# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:42:49 2018

This module includes classes and functions relating arrays.

@author: Guest Group
"""
from warnings import warn
from biosteam import np
from copy import copy

__all__ = ('array', 'tuple_array', 'material_array', 'mass_array', 'organized_list', 'same_type_tuple', 'same_type', 'reorder', 'get_frac')

ndarray = np.ndarray
asarray = np.asarray
integer = np.integer

# %% New Array Like Classes

class array(ndarray):
    """Create a numpy ndarray."""
    # This allows casting to a module array class
    __slots__ = ()

class tuple_array(ndarray):
    """Create an array that is immutable and hashable.

    **Parameters**

         **array:** [array_like]
         
         **dtype:** [numpy datatype]
         
         **order:** {'C', 'F'}

    """
    __slots__ = ()
    
    def __new__(cls, arr, dtype=np.float64, order=None):
        self = asarray(arr, dtype, order).view(cls)
        return self

    def __hash__(self):
        return hash((ndarray,) + tuple(self))

    def __setitem__(self, key, val):
        raise TypeError(f"'{type(self).__name__}' objects are immutable.")

    def __array_wrap__(self, out_arr):
        # New array should not be a tuple_array
        if self is out_arr:
            raise TypeError(f"'{type(self).__name__}' objects are immutable.")
        out_arr.__class__ = array
        return out_arr


class material_array(ndarray):
    """Create an array that issues a RuntimeWarning when a non-positive or non-finite value is encountered.

    **Parameters**

         **array:** [array_like]
         
         **dtype:** [numpy datatype]
         
         **order:** {'C', 'F'}

    """
    __slots__ = ()
    @classmethod
    def enforce_valuecheck(cls, val):
        """If *val* is True, issue warning when non-positive or non-finite values are encountered."""
        if val:
            cls.__setitem__ = cls._setitem 
            cls.__array_wrap__ = cls._array_wrap
        else:
            cls.__setitem__ = array.__setitem__
            cls.__array_wrap__ = array.__array_wrap__

    def __new__(cls, arr, dtype=np.float64, order=None):
        self = asarray(arr, dtype, order).view(cls)
        return self

    def __setitem__(self, key, val):
        # When self[:] = self
        if val is self: return
        
        # Check values and return
        val = np.asarray(val)
        if (val < 0).any() or (val == np.inf).any():
            invalid = val[(val <= 0) | (val == np.inf)]
            warn(RuntimeWarning(f"Encountered non-positive or non-finite value {invalid} in '{type(self).__name__}' object."), stacklevel=2)
        super().__setitem__(key, val)
    
    _setitem = __setitem__
    
    def __array_wrap__(self, out_arr):
        # New array should not be a material_array
        if self is out_arr:
            if (self < 0).any() or (self == np.inf).any():
                invalid = self[(self < 0) | (self == np.inf)]
                warn(RuntimeWarning(f"Encountered non-positive or non-finite value {invalid} in '{type(self).__name__}' object."), stacklevel=2)
        else:
            out_arr.__class__ = array
        return out_arr
    
    _array_wrap = __array_wrap__


class mass_array(ndarray):
    __slots__ = ('MW', 'array', '_array')
    
    def __new__(cls, MW, array):
        self = asarray(MW*array, dtype=np.float64).view(cls)
        self.MW = MW
        self.array = array
        self._array = copy(array)
        return self
    
    def __getitem__(self, key):
        MW = self.MW[key]
        last = self._array[key]
        current = self.array[key]
        if isinstance(key, (int, integer)):
            if last == current:
                value = super().__getitem__(key)
            else:
                value = current*MW
                super().__setitem__(key, value)
        else:
            value = super().__getitem__(key)
            i = 0
            for L, C  in zip(last, current):
                if L != C:
                    super(mass_array, value).__setitem__(i, C*MW[i])
                    last[i] = C
                i += 1            
            value.MW = MW
            value.array = current
            value._array = last
        return value
            
    def __setitem__(self, key, value):
        MW = self.MW[key]
        slice_self = super().__getitem__(key)
        if isinstance(key, (int, integer)):
            if value != slice_self:
                self.array[key] = self._array[key] = value/MW
                super().__setitem__(key, value)
            return
        last = self.array[key]
        current = self._array[key]
        i = 0
        for v in value:
            a = super(mass_array, slice_self).__getitem__(i)
            if v != a:
                last[i] = current[i] = v/MW[i]
                super(mass_array, slice_self).__setitem__(i, v)
            i += 1
    
    def __array_wrap__(self, out_arr):
        # Unrelated array should not be a mass_array
        if not hasattr(out_arr, 'MW'):
            out_arr.__class__ = not_mass_array
        return out_arr
    
    def __str__(self):
        return str(np.asarray(self[:]))
        
    def __repr__(self):
        return 'mass_' + repr(np.asarray(self[:]))

class not_mass_array(ndarray):
    __slots__ = mass_array.__slots__
    
not_mass_array.__name__ = 'array'

# %% Functions for array like objects

def reorder(array: 'array_like', current_order: tuple, new_order: tuple) -> np.ndarray:
    """Return a reordered array with zero in place of elements that exist in the new order but not in the current order.
    
    .. code-block:: python
       
       >>> reorder([1,2,3], ('a', 'b', 'c'), ('a', 'c', 'd', 'b'))
       array([1, 3, 0, 2])
    
    """
    n = 0
    p_new = np.zeros(len(new_order))
    for ID1 in current_order:
        m = 0
        for ID2 in new_order:
            if ID1 == ID2:
                p_new[m] = array[n]
            m += 1
        n += 1
    return p_new

def get_frac(array) -> 'Normalized array':
    """Normalize array to a magnitude of 1. If magnitude is zero, all fractions will have equal value."""
    sum_array = sum(array)
    if sum_array == 0:
        l = len(array)
        frac = np.ones(l)/l
    else:
        frac = array/sum_array
    return frac

def organized_list(iterable):
    """Return a list with iterable in alphabetical order using its string representation. Repetitions are not included."""
    iterable = sorted(set(i for i in iterable if str(i) != ''), key=lambda i: str(i))
    return iterable

def same_type_tuple(type_, iterable):
    """Return the iterable as a tuple. Raise TypeError if any item in iterable is not an instance of type_."""
    # Make sure iterable are in a tuple
    if isinstance(iterable, type_):
        iterable = (iterable,)
    else:
        iterable = tuple(iterable)    
        # Check that all are type_ instances
        for s in iterable:
            if not isinstance(s, type_):
                raise TypeError(f"Only '{type_.__name__}' objects are valid elements for ins, not '{type(s).__name__}' objects.")
    return iterable

def same_type(iterable, type_):
    """Raise TypeError if any item in iterable is not an instance of type_."""
    # Make sure iterable are in a tuple
    if not isinstance(iterable, type_):
        # Check that all are type_ instances
        for s in iterable:
            if not isinstance(s, type_):
                raise TypeError(f"Only '{type_.__name__}' objects are valid elements, not '{type(s).__name__}' objects.")