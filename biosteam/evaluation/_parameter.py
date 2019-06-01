# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:20:53 2019

@author: yoelr
"""
__all__ = ('Parameter',)
from ._name import elementname

class Parameter:
    __slots__ = ('name', 'setter', 'simulate', 'element', 'system')
    
    def __init__(self, name, setter, simulate, element, system):
        self.name = name.replace('_', ' ').capitalize()
        self.setter = setter
        self.simulate = simulate
        self.element = element
        self.system = system
    
    @property
    def element_name(self):
        return elementname(self.element)
    
    def __repr__(self):
        if self.element:
            return f'<{type(self).__name__}: [{self.element_name}] {self.name}>'
        else:
            return f'<{type(self).__name__}: {self.name}>'
    
               