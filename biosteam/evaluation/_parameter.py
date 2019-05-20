# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:20:53 2019

@author: yoelr
"""
__all__ = ('Parameter',)

class Parameter:
    __slots__ = ('_name', '_setter', '_simulate', '_element', '_system')
    
    def __init__(self, name, setter, simulate, element, system):
        self._name = name.replace('_', ' ').capitalize()
        self._setter = setter
        self._simulate = simulate
        self._element = element
        self._system = system
    
    def update(self, value):
        self._setter(value)
        self._simulate()
    
    @property
    def name(self):
        return self._name
    
    def __repr__(self):
        return f'<{type(self).__name__}: {self._name}>'
    
               