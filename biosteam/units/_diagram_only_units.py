# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 23:25:42 2020

@author: yoelr
"""
from biosteam._unit import Unit
from .._graphics import system_unit, stream_unit

__all__ = ('DiagramOnlyUnit',
           'DiagramOnlySystemUnit',
           'DiagramOnlyStreamUnit')

class DiagramOnlyUnit(Unit, isabstract=True):
    _ID = ID = None
    _N_ins = _N_outs = 1
    _ins_size_is_fixed = _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._register(ID)
    
    def _register(self, ID): 
        self.ID = self._ID = ID
    
class DiagramOnlySystemUnit(DiagramOnlyUnit, isabstract=True):
    """Dummy unit for displaying a system as a unit."""
    line = 'System'
    _graphics = system_unit

class DiagramOnlyStreamUnit(DiagramOnlyUnit, isabstract=True):
    """Dummy unit for displaying a streams as a unit."""
    line = ''
    _graphics = stream_unit