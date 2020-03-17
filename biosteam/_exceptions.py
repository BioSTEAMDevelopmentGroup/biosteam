#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:49:52 2018

This module includes classes and functions relating exception handling.

@author: Yoel Rene Cortes-Pena
"""
from .utils.biosteam_colors import colors

__all__ = ('DesignError',)

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""

def message_with_object_stamp(object, msg):
    return colors.violet(repr(object)) + ' ' + msg

def raise_error_with_object_stamp(object, error):
    if hasattr(error, 'args'):
        msg, *args = error.args
        error.args = (message_with_object_stamp(object, msg), *args)
    if hasattr(error, 'msg'):
        error.msg = message_with_object_stamp(object, error.msg)
    raise error

def try_method_with_object_stamp(object, method):
    try:
        method()
    except Exception as error:
        raise_error_with_object_stamp(object, error)
