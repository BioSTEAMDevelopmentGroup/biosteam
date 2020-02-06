# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 07:19:49 2019

@author: yoelr
"""
from warnings import warn

__all__ = ('DesignWarning', 'design_warning',
           'lb_warning', 'ub_warning', 'bounds_warning')

class DesignWarning(Warning):
    """Warning regarding design constraints."""

# %% Bounds checking

def design_warning(source, msg):
    """Return a DesignWarning object with source description."""
    if isinstance(source, str):
        msg = f'@{source}: ' + msg
    elif source:
        msg = f'@{type(source).__name__} {str(source)}: ' + msg
    return DesignWarning(msg)
            
def lb_warning(key, value, units, lb, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb}{units})."
    
    warn(design_warning(source, msg, DesignWarning), stacklevel=stacklevel)
    
def ub_warning(key, value, units, ub, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub}{units})."
    
    warn(design_warning(source, msg, DesignWarning), stacklevel=stacklevel)

def bounds_warning(source, key, value, units, bounds):
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
        warn(design_warning(source, msg), stacklevel=3)