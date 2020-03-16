# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 02:24:35 2019

@author: yoelr
"""
from ._unit import Unit

__all__ = ('Facility',)

def get_network_priority(facility):
    return facility.network_priority

class Facility(Unit, isabstract=True):
    
    @staticmethod
    def ordered_facilities(facilities):
        """Return facilitied ordered according to their network priority."""
        return sorted(facilities, key=get_network_priority)
    
    def __init_subclass__(cls, isabstract=False):
        super().__init_subclass__(isabstract)
        if not hasattr(cls, 'network_priority'):
            raise NotImplementedError('Facility subclasses must implement a '
                                      '`network_priority` attribute to designate '
                                      'the order of simulation relative to other '
                                      'facilities')
    
    @property
    def system(self):
        return self._system
        
