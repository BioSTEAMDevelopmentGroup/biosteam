# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from ._unit import Unit
from .exceptions import UnitInheritanceError

__all__ = ('Facility',)

def get_network_priority(facility):
    return facility.network_priority

class Facility(Unit, isabstract=True,
               new_graphics=False):
    """Abstract class for facilities that are run after simulation of all
    unit operations within a system path."""
    _universal = True
    autonumber = False # Default ID will not include number
    _skip_simulation_when_inlets_are_empty = False # Should be false with most facilities.
    @staticmethod
    def ordered_facilities(facilities):
        """Return facilities ordered according to their network priority."""
        return sorted(set(facilities), key=get_network_priority)
    
    def __init_subclass__(cls,
                          isabstract=False,
                          new_graphics=True):
        super().__init_subclass__(isabstract,
                                  new_graphics)
        if not hasattr(cls, 'network_priority'):
            raise UnitInheritanceError(
                'Facility subclasses must implement a `network_priority` '
                'attribute to designate the order of simulation relative to '
                'other facilities'
            )
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, **kwargs):
        Unit.__init__(self, ID, ins, outs, thermo, **kwargs)
        self._system = None
        self._other_units = None
    
    @property
    def other_units(self):
        return self._other_units
        