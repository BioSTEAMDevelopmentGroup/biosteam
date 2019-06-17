# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 00:53:29 2019

@author: yoelr
"""
from ... import Unit

__all__ = ('extend_finalize', )

def extend_finalize(cls):
    """Extends the Unit class with the following abstract methods:

    **_spec()**
        Calculate specification cost factors and update purchase prices (called after `_cost`).
        
    **_end():**
        Finish setting purchase prices and utility costs (called after `_spec`).
        
    """

    if cls._finalize is Unit._finalize:
        cls._finalize = _finalize
    elif cls._finalize is not _finalize:
        raise RuntimeError("cannot decorate Unit class with an implemented '_finalize' method")

def _finalize(self):
    """Run all cost methods and finalize costs."""
    self._cost()
    self._spec()
    self._end()
    self._update_capital_cost()
    self._update_utility_cost()