# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:42:49 2018

This module includes classes and functions relating arrays.

@author: Guest Group
"""
import numpy as np

__all__ = ('fraction',)

ndarray = np.ndarray
asarray = np.asarray

# %% Functions

def fraction(array):
    """Return a normalized array to a magnitude of 1. If magnitude is zero, all fractions will have equal value."""
    sum_array = array.sum()
    if sum_array == 0:
        l = len(array)
        return np.ones(l)/l
    else:
        return array/sum_array
    