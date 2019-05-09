# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:22:37 2019

@author: yoelr
"""
        
from .. import Unit, Stream
from inspect import signature
from .._utils import function

__all__ = ('Block',)

do_nothing = lambda: None
    

# %% Simulation blocks

class Block:
    """
    Create a Block object that can simulate the element and the system downstream. The block can also generate block functions.
    
    **Parameters**
    
        **element:** [Unit or Stream] Element in the system.
        
        **system:** [System] If None, ignore downstream.
        
    **Examples**

         :doc:`Block Example`
    
    Create a block object passing a Unit object, and simulate the block:
        
    .. code-block:: python
    
       >>> P1 = Unit('P1')
       >>> block = Block(P1)
       >>> block.show()
       [Block: Unit-P1]
       >>> block.simulate()
    
    .. Note::
    
       Because no system was passed, the block would simulate just the unit.
       
    Create a block object passing both a Unit object and a System object:
        
    .. code-block:: python
    
       >>> P0 = Unit('P0', outs=Stream())
       >>> P1 = Unit('P1', ins=P0-0)
       >>> P2 = Unit('P2', ins=P1-0)
       >>> system = System('', network=(P0, P1, P2))
       >>> block = Block(P1, system)
       >>> block.show()
       [Block: Unit-P1 and downstream]
       >>> block.system.show()
       System: Unit-P1 and downstream
        network: (P1, P2)
    
    .. Note::
    
       The block would simulate Unit-P1 and downstream.
    
    
    Calling a block object with a setter will return a block function equivalent to:
        
    .. code-block:: python
    
       >>> def blockfunc(args):
       ...     setter(args)
       ...     self.simulate() # self is the Block object
       
    For example:
        
    .. code-block:: python
    
       >>> # setter functions is hypothetical
       >>> def setter(args): pass
       >>> blockfunc = Block(P1, None)(setter)
       >>> blockfunc
       <function [Block P1] setter(args)>
    
    .. Note::
        
       The function name and signature matches the setter function.
    
    """
    
    _cachedblocks = {}
    __slots__ = ('_system', '_simulate', '_element')
    
    def __new__(cls, element, system=None):
        cached = cls._cachedblocks
        block = cached.get((system, element))
        if block:
            self = block
        else:
            self = object.__new__(cls)
            self.__init__(element, system)
        return self
    
    def __init__(self, element, system=None):
        inst = isinstance
        if inst(element, Stream): unit = element.sink
        elif inst(element, Unit): unit = element
        if system and element:
            subsys = system._downstream_system(unit)
            simulate = subsys.simulate
        else:
            subsys = system
            simulate = unit.simulate if inst(element, Unit) else do_nothing
            
        self._system = subsys
        self._simulate = simulate
        self._element = element
        self._cachedblocks[system, element] = self
    
    def __call__(self, setter, simulate=None, param=None) -> function:
        """Return a block function."""
        if simulate is None: simulate = self._simulate
        if not param: param, = signature(setter).parameters.keys()
        
        # Make blockfunc
        name = setter.__name__
        if name[0] == '<': name = 'Lambda'
        str2exec =  (f'def {name}({param}):\n'
                   + f'    setter({param})\n'
                   + f'    simulate()')

        globs = {'setter': setter,
                 'simulate': simulate}
        locs = {}
        exec(str2exec, globs, locs)
        blockfunc = locs[name]
        blockfunc.__qualname__ = f'{self} {blockfunc.__qualname__}'
        blockfunc._param = param
        blockfunc._simulate = simulate
        blockfunc._element = self._element
        blockfunc._system = self._system
        blockfunc._setter = setter
        return blockfunc
    
    @property
    def system(self):
        """System that the block simulates."""
        return self._system
    
    @property
    def element(self):
        """Block element that gets checked for changes before simulation."""
        return self._element
    
    @property
    def simulate(self):
        """Simulate block system if block element has changed. If block element has not changed, simulate only the element."""
        return self._simulate
    
    #: A Block object without a Unit or a System object.
    _emptyblock = None
        
    def _repr(self):
        if self._element:
            if isinstance(self._element, type):
                return f"{type(self).__name__}: {self._element.__name__}"
            return f"{type(self).__name__}: {type(self._element).__name__}-{self._element}"
        else:
            return f"{type(self).__name__}"
    
    def __repr__(self):
        return f'[{self._repr()}]'
    
    def _info(self):
        if self._system:
            return f'[{type(self).__name__}: {self._system}]'
        else:
            return repr(self)

    def show(self):
        print(self._info())
