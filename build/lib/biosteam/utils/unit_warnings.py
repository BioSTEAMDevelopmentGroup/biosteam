# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 07:19:49 2019

@author: yoelr
"""
from warnings import warn
from .._exceptions import message_with_object_stamp

__all__ = ('UnitWarning',
           'DesignWarning', 
           'CostWarning',
           'lb_warning',
           'ub_warning',
           'bounds_warning')

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

def design_warning_with_source(source, msg):
    """Return a DesignWarning object with source description."""
    msg= message_with_object_stamp(source, msg)
    return DesignWarning(msg)
            
def lb_warning(key, value, units, lb, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb}{units})."
    
    warn(DesignWarning.from_source(source, msg), stacklevel=stacklevel)
    
def ub_warning(key, value, units, ub, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub}{units})."
    
    warn(DesignWarning.from_source(source, msg), stacklevel=stacklevel)

def bounds_warning(source, key, value, units, bounds, kind='design'):
    """Issue a warning if value is out of bounds.
    
    Parameters
    ----------
    key : str
          Name of value.
    value : float
    units : str
            Units of value        
    bounds : iterable[float, float]
             Upper and lower bounds.
        
    """
    # Warn when value is out of bounds
    lb, ub = bounds
    if not (lb <= value <= ub):
        units = ' ' + units if units else ''
        try:
            msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb:.4g} to {ub:.4g}{units})."
        except:  # Handle format errors
            msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb} to {ub}{units})."
        if kind == 'design':
            Warning = DesignWarning
        elif kind == 'cost':
            Warning = CostWarning
        else:
            raise ValueError(f"kind must be either 'design' or 'cost', not {repr(kind)}")
        warning = Warning.from_source(source, msg)
        warn(warning, stacklevel=3)