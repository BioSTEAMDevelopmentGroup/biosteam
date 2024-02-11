# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import exceptions
from thermosteam.exceptions import *
from warnings import warn

__all__ = (
    'DesignError', 
    'FailedEvaluation',
    'Converged',
    'UnitWarning',
    'DesignWarning', 
    'CostWarning',
    'lb_warning',
    'ub_warning',
    'bounds_warning',
    *exceptions.__all__)
del exceptions

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""
    
class FailedEvaluation(RuntimeWarning):
    """RuntimeWarning regarding failed model evaluation."""
    
class Converged(Exception):
    """Exception to stop iteration early."""
    
# %% BioSTEAM warnings

class UnitWarning(Warning):
    """Warning regarding unit operations."""
    
    @classmethod
    def from_source(cls, source, msg):
        """Return a DesignWarning object with source description."""
        msg= message_with_object_stamp(source, msg)
        return cls(msg)
    
class DesignWarning(UnitWarning):
    """Warning regarding design constraints."""
    
class CostWarning(UnitWarning):
    """Warning regarding design constraints."""

# %% Bounds checking

def design_warning_with_source(source, msg): # pragma: no cover
    """Return a DesignWarning object with source description."""
    msg= message_with_object_stamp(source, msg)
    return DesignWarning(msg)
            
def lb_warning(source, key, value, units, lb, stacklevel=2): # pragma: no cover
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb}{units})."
    
    warn(DesignWarning.from_source(source, msg), stacklevel=stacklevel)
    
def ub_warning(source, key, value, units, ub, stacklevel=2): # pragma: no cover
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub}{units})."
    
    warn(DesignWarning.from_source(source, msg), stacklevel=stacklevel)

def bounds_warning(source, key, value, units, bounds, kind='design'): # pragma: no cover
    """Issue a warning if value is out of bounds.
    
    Parameters
    ----------
    source : Unit
        Unit where the warning is issued
    key : str
        Name of value.
    value : float
    units : str
        Units of value        
    bounds : Iterable[float, float]
        Upper and lower bounds.
        
    """
    # Warn when value is out of bounds
    lb, ub = bounds
    if not (lb <= value <= ub):
        units = ' ' + units if units else ''
        try:
            msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb:.4g} to {ub:.4g}{units})"
        except:  # Handle format errors
            msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb} to {ub}{units})"
        if kind == 'design':
            Warning = DesignWarning
            msg += ' for design algorithm'
        elif kind == 'cost':
            msg += ' for cost correlation'
            Warning = CostWarning
        else:
            raise ValueError(f"kind must be either 'design' or 'cost', not {repr(kind)}")
        warning = Warning.from_source(source, msg)
        warn(warning, stacklevel=3)