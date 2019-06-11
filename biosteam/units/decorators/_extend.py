# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 00:53:29 2019

extend_summary extends the Unit with the following abstract methods:

**_N()**
    Add number of units to "Design" dictionary (called after `_design`).

extend_finalize extends the Unit with the following abstract methods:

**_spec()**
    Calculate specification cost factors and update purchase prices (called after `_cost`).
    
**_end():**
    Finish setting purchase prices and utility costs (called after `_spec`).


@author: yoelr
"""
from ... import Unit

__all__ = ('extend_finalize', 'extend_summary')

def extend_summary(cls):
    if cls._summary is Unit._summary:
            cls._summary = _summary
    elif cls._summary is not _summary:
        raise RuntimeError("cannot decorate Unit class with an implemented '_summary' method")

def extend_finalize(cls):
    if cls._finalize is Unit._finalize:
        cls._finalize = _finalize
    elif cls._finalize is not _finalize:
        raise RuntimeError("cannot decorate Unit class with an implemented '_summary' method")

def _summary(self):
    """Calculate all results from unit run."""
    self._design()
    self._N()
    self._finalize()

def _finalize(self):
    """Run all cost methods and finalize purchase and utility cost."""
    self._cost()
    self._spec()
    self._end()
    self._update_utility_cost()
    self._update_purchase_cost()