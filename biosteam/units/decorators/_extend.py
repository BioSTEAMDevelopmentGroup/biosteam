# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 00:53:29 2019

@author: yoelr
"""
from ... import Unit

__all__ = ('extend_summary', )

def extend_summary(cls):
    """Extends the Unit class with the following abstract methods:
        
    **_end():**
        Finish setting purchase prices and utility costs.
        
    """
    if hasattr(cls, '_end'):    
        if cls._summary is Unit._summary:
            cls._summary = _summary
        elif cls._summary is not _summary:
            raise RuntimeError("cannot decorate Unit subclass an implemented '_summary' method")

def _summary(self):
    """Calculate all results from unit run."""
    self._design()
    self._cost()
    self._end()