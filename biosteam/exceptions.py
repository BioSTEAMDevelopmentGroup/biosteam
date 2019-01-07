#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:49:52 2018

This module includes classes and functions relating exception handling.

@author: Yoel Rene Cortes-Pena
"""

from warnings import warn

# %% Biosteam errors

class biosteam_Error(Exception):
    """General exception in biosteam."""


class Stop(biosteam_Error):
    """Exception to Stop code."""


class IDconflict(biosteam_Error):
    """Exception regarding ID conflicts."""


class DesignError(biosteam_Error):
    """Exception regarding unit design."""


class metaError(biosteam_Error):
    """Exception regarding metaclasses."""


class SolverError(biosteam_Error):
    """Exception regarding solvers."""


class EquilibriumError(biosteam_Error):
    """Exception regarding equilibrium."""


class DimensionError(biosteam_Error):
    """Exception regarding wrong dimensions."""


# Python's built in KeyError quotes the message, used as a by-pass for the debbuger
KE = type('KeyError', (biosteam_Error, ), {})


# %% Exception handling

def ignore_wrapper(function):
    """Wrap function with a try-except clause to ignore any errors."""
    def ignore_function(*args, **kwargs):
        try:
            function(*args, **kwargs)
        except:
            pass
    return ignore_function



