# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import Unit
from thermosteam import Stream
from inspect import signature
from ._parameter import Parameter

__all__ = ('Block',)

def do_nothing(): pass

# %% Simulation blocks

class Block:
    """
    Create a Block object that can simulate the element and the system 
    downstream. The block can also generate Parameter objects that can
    update the system state.
    
    Parameters
    ----------
    element : Unit or :class:`~thermosteam.Stream`
        Element in the system.  
    system=None : System, optional
        System affected by element. In other words, element should affect
        downstream operations in the system.
        
    """
    
    _blocks = {}
    __slots__ = ('_system', '_simulate', '_element')
    
    def __new__(cls, element, system=None):
        block = cls._blocks.get((system, element))
        if block:
            self = block
        else:
            self = object.__new__(cls)
            self.__init__(element, system)
        return self
    
    def __init__(self, element, system=None):
        isa = isinstance
        if isa(element, Stream): unit = element.sink
        elif isa(element, Unit): unit = element
        if system and element:
            subsys = system._downstream_system(unit)
            simulate = subsys.simulate
        else:
            subsys = system
            simulate = unit.simulate if isa(element, Unit) else do_nothing
        self._system = subsys
        self._simulate = simulate
        self._element = element
        self._blocks[system, element] = self
    
    def parameter(self, setter, simulate, name, distribution, units, 
                  baseline, bounds) -> Parameter:
        """Return a Parameter object."""
        if simulate is None: simulate = self._simulate
        if not name: name, = signature(setter).parameters.keys()
        return Parameter(name, setter, simulate, self._element,
                         self._system, distribution, units, baseline, bounds)
    
    @property
    def system(self):
        """System that the block simulates."""
        return self._system
    
    @property
    def element(self):
        """Starting element of the system (if any)."""
        return self._element
    
    @property
    def simulate(self):
        """Simulate block."""
        return self._simulate
        
    def __repr__(self):
        if self._system:
            return f'<{type(self).__name__}: {self._system}>'
        elif self._element:
            if isinstance(self._element, type):
                return f"<{type(self).__name__}: {self._element.__name__}>"
            return f"<{type(self).__name__}: {type(self._element).__name__}-{self._element}>"
        else:
            return f"<{type(self).__name__}>"
    
    
    
    
    