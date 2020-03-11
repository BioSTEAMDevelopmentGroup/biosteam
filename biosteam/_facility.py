# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from ._unit import Unit

__all__ = ('Facility',)

class Facility(Unit, isabstract=True, new_graphics=False):
    
    _is_facility_ = True
    
    @property
    def system(self):
        return self._system
        