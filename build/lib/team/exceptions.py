#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:49:52 2018

This module includes classes and functions relating exception handling.

@author: Yoel Rene Cortes-Pena
"""
import sys

# %% Biosteam errors

class biosteamError(Exception):
    """General exception in biosteam."""

class Stop(biosteamError):
    """Exception to Stop code."""

class DesignError(biosteamError):
    """Exception regarding unit design."""

class SolverError(biosteamError):
    """Exception regarding solvers."""

class EquilibriumError(biosteamError):
    """Exception regarding equilibrium."""

class DimensionError(biosteamError):
    """Exception regarding wrong dimensions."""


# Python's built in KeyError quotes the message, used as a by-pass for the debbuger
KE = type('KeyError', (biosteamError, ), {})


# %% Biosteam Warnings

class DesignWarning(Warning):
    """Warning regarding design constraints."""
    
#%% Decorators and functions

def notify_error(func):
    """Decorate class method to provide a location summary when an error occurs."""
    if hasattr(func, '_original'):
        func = func._original

    def wrapper(self, *args, **kwargs):
        try: return func(self, *args, **kwargs)
        except Exception as e:
            # If exception already include location, it is replaced
            location = f'@{type(self).__name__} {self}'
            msg = str(e).strip('\n').replace(location + ': ', '')
            
            # Add location to message
            if not ('@' in msg and ':\n' in msg):
                msg = location + f'.{func.__name__}:\n' + msg                 
            
            # Raise exception with same traceback but new message
            if type(e) is KeyError:
                raise KE(msg).with_traceback(sys.exc_info()[2])
            else:
                raise type(e)(msg).with_traceback(sys.exc_info()[2])
    
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    wrapper.__annotations__ = func.__annotations__
    return wrapper



