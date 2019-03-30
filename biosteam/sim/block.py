# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:22:37 2019

@author: yoelr
"""
        
from biosteam import Unit, Stream, np
from inspect import signature
from biosteam.utils import function
from copy import copy

__all__ = ('Block', 'emptyblock')

do_nothing = lambda: None
    

# %% Simulation blocks

class Block:
    """
    Create a Block object that can simulate the element and the system downstream. The block can also generate block functions.
    
    **Parameters**
    
        **element:** [Unit or Stream] Element in the system.
        
        **system:** [System] If None, ignore downstream.
        
        **getter:** [function] Getter for block function.
        
        **setter:** [function] Setter for block function.
        
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
    
    
    Calling a block object with a getter and a setter will return a block function equivalent to:
        
    .. code-block:: python
    
       >>> def blockfunc(args):
       ...     setter(args)
       ...     self.simulate() # self is the Block object
       ...     return getter()  
       
    For example:
        
    .. code-block:: python
    
       >>> # getter and setter functions are hypothetical
       >>> def getter(): pass
       >>> def setter(args): pass
       >>> blockfunc = Block(P1, None)(getter, setter)
       >>> blockfunc
       <function [Block P1] getter(args)>
    
    .. Note::
        
       The function signature matches the setter while the function name matches the getter.
    
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
        else: raise ValueError(f"element must be either a Unit, a Stream object, not '{type(element).__name__}'.")
        if system:
            subsys = system._downstream_system(unit)
            simulate = subsys.simulate
        else:
            subsys = system
            if element is Unit: simulate = unit.simulate
            else: simulate = do_nothing
            
        self._system = subsys
        self._simulate = simulate
        self._element = element
        self._cachedblocks[system, element] = self
    
    def _make_blockfunc(self, getter, setter, docfunc):
        # Prevent downstream exceptions
        getter_param = signature(getter).parameters
        if getter_param:
            for p in getter_param.values():
                if '=' not in str(p):
                    raise ValueError(f'{type(self).__name__} getter signature cannot have arguments.')
        sig = signature(setter)
        setter_params = tuple(sig.parameters.keys())
        param = ', '.join(setter_params)
        if len(setter_params) != 1: 
            raise ValueError(f'{type(self).__name__} setter signature must have one and only one argument.')
        
        # Make blockfunc
        name = getter.__name__
        if name[0] == '<': name = 'Lambda'
        str2exec =  (f'def {name}{sig}:    \n'
                   + f'    setter({param}) \n'
                   + f'    simulate()      \n'
                   + f'    return getter()   ')

        globs = {'setter': setter,
                 'simulate': self._simulate,
                 'getter': getter}
        locs = {}
        exec(str2exec, globs, locs)
        blockfunc = locs[name]
        blockfunc.__qualname__ = f'{self} {blockfunc.__qualname__}'
        blockfunc.__doc__ = docfunc.__doc__
        blockfunc._param = param
        blockfunc._block = self
        blockfunc._setter = setter
        blockfunc._getter = getter
        return blockfunc
    
    def __call__(self, getter, setter) -> function:
        return self._make_blockfunc(getter, setter, getter)
    
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
    emptyblock = None
        
    def _repr(self):
        if self._element:
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

emptyblock = object.__new__(Block)
emptyblock._system = None
emptyblock._element = None
emptyblock._simulate = do_nothing
Block.emptyblock = emptyblock