#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:49:52 2018

This module includes classes and functions relating exception handling.

@author: Yoel Rene Cortes-Pena
"""
# import sys as _sys

__all__ = ('DesignError',)

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""

# class UnitError(RuntimeError):
#     """RuntimeError regarding unit operations."""
#     def __init__(self, unit, method, error):
#         # Add location to message
#         msg = (f'{type(error).__name__} at {repr(unit)}.{method.__name__}\n{error}')
        
#         # Raise exception with same traceback but new message
#         super().__init__(msg)
#         self.original_exception = error
#         self.with_traceback(_sys.exc_info()[2])


#%% Decorators and functions

# def _try_method(method):
#     try: return method()
#     except UnitError:
#         raise UnitError
#     except Exception as error:
#         # If exception already include location, it is replaced
#         try:
#             self = method.__self__
#         except:
#             raise error
#         raise UnitError(self, method, error)
        

