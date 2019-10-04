# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from ._unit import Unit
from .utils import NotImplementedMethod

__all__ = ('Facility',)

class Facility(Unit):
    
    _run = NotImplementedMethod
    
    @property
    def system(self):
        return self._system
        