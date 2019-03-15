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

    

# %% Simulation blocks

class Block:
    """Create a Block object that can simulate an `element` and the `system` downstream.
    
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
    
       Simulating the block will simulate the unit only if kwargs of the unit has changed. Otherwise, it will only cost the unit.
       
    Create a block object passing both a Unit object and a System object:
        
    .. code-block:: python
    
       >>> P0 = Unit('P0', outs=Stream())
       >>> P1 = Unit('P1', ins=P0-0)
       >>> P2 = Unit('P2', ins=P1-0)
       >>> system = System('', network=[P0, P1, P2])
       >>> block = Block(P1, system)
       >>> block.show()
       [Block: Unit-P1 and downstream]
       >>> block.system.show()
       System: Unit-P1 and downstream
        network: (P1, P2)
       
    .. Note::
        
       Simulating the block will simulate the system only if kwargs of the unit has changed. Otherwise, it wil only cost the system.
    
    If both a getter and a setter are given, create a block function equivalent to:
        
    .. code-block:: python
    
       >>> def blockfunction(args):
       ...     setter(args)
       ...     self.simulate() # self is the Block object
       ...     return getter()  
       
    For example:
        
    .. code-block:: python
    
       >>> def getter(): pass
       >>> def setter(efficiency): P1.kwargs['nu'] = efficiency
       >>> block_function = Block(P1, None, getter, setter)
       >>> block_function
       <function [Block P1] getter(efficiency)>
    
    .. Note::
        
       The function signature matches the setter while the function name matches the getter.
    
    If either a getter or a setter is given, create a decorator that returns a block function with the missing argument. The overall format is:
    
    .. code-block:: python
    
       >>> @Block(P1, system, getter)
       ... def func(args): pass
       >>> func
       <function [Block P1] getter(args)>
       
    Alternatively:
        
    .. code-block:: python
    
       >>> @Block(P1, system, setter=setter)
       ... def func(args): pass
       >>> func
       <function [Block P1] func(args)>
    
    Block objects serve as decorators that return a block function given a getter or a setter:
        
    .. code-block:: python
       
       >>> block = Block(P1, system)
       >>> @block(getter)
       ... def func(args): pass
       >>> func
       <function [Block P1] getter(args)>
    
    """
    
    _cachedblocks = {}
    __slots__ = ('_system', '_simulate', '_element')
    
    def __new__(cls, element, system=None, getter=None, setter=None):
        cached = cls._cachedblocks
        block = cached.get((system, element))
        if block:
            self = block
        else:
            self = object.__new__(cls)
            self.__init__(element, system)
        if getter or setter: return self(getter, setter)
        else: return self
    
    def __init__(self, element, system=None):
        if isinstance(element, Stream):
            unit = element.sink[0]
        elif isinstance(element, Unit):
            unit = element
        else:
            raise ValueError(f"element must be either a Unit, a Stream object, not '{type(element).__name__}'.")
        if not system:
            subsys = system
            if element is Unit:
                def simulate(): unit.simulate()
            else:
                def simulate(): pass
        
        elif element is unit:
            subsys = system._downstream_system(unit)
            def simulate():
                if element._kwargs == element.kwargs:
                    element._summary()
                else:
                    element._kwargs = copy(element.kwargs)
                    subsys._reset_iter()
                    unit._setup()
                    subsys._converge()
                    for u in subsys.units:
                        u._summary()
                for i in subsys.facilities:
                    if isinstance(i, function):
                        i()
                    else:
                        i._run()
                        i._summary()
        else: # element is a stream
            subsys = system._downstream_system(unit)
            molarray = element._molarray
            molarray_old = np.array(molarray)
            T_old = element.T
            P_old = element.P
            def simulate():
                nonlocal molarray_old, T_old, P_old
                if ((molarray == molarray_old).all()
                    and element.T == T_old
                    and element.P == P_old):
                    return
                subsys._reset_iter()
                subsys._converge()
                for u in subsys.units:
                    u._summary()
                for i in subsys.facilities:
                    if isinstance(i, function):
                        i()
                    else:
                        i._run()
                        i._summary()
                molarray_old = np.array(molarray)
                T_old = element.T
                P_old = element.P
        self._system = subsys
        self._simulate = simulate
        self._element = element
        self._cachedblocks[system, element] = self
    
    def _make_blockfunc(self, getter, setter, docfunc):
        # Prevent downstream exceptions
        getter_param =signature(getter).parameters
        if getter_param:
            for p in getter_param.values():
                if '=' not in str(p):
                    raise ValueError(f'{type(self).__name__} getter signature cannot have arguments.')
        sig = signature(setter)
        # import pdb
        # pdb.set_trace()
        setter_params = tuple(sig.parameters.keys())
        param = ', '.join(setter_params)
        if len(setter_params) != 1: 
            raise ValueError(f'{type(self).__name__} setter signature must have one and only one argument.')
        
        # Make blockfunc
        name = getter.__name__
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
        # annotations = blockfunc.__annotations__
        # if 'return' in getter.__annotations__:
        #     annotations['return'] = getter.__annotations__['results']
        # for i in setter_params:
        #     sa = setter.__annotations__
        #     if i in sa: annotations[i] = sa[i]
        return blockfunc
    
    def _make_decorator(self, getter, setter):
        globs = {'make_blockfunc': self._make_blockfunc,
                 'getter': getter,
                 'setter': setter,}
        sig = 'setter' if getter else 'getter'
        str2exec =  (
              f'def decorator({sig}):                                   \n'
            + f'    """Decorate as a {type(self).__name__} function.""" \n'
            + f'    return make_blockfunc(getter, setter, {sig})          ')
        locs = {}
        exec(str2exec, globs, locs)
        decorator = locs['decorator']
        decorator.__qualname__ = f'{self} {decorator.__qualname__}'
        return decorator
    
    def __call__(self, getter=None, setter=None) -> function:
        if getter and setter:
            return self._make_blockfunc(getter, setter, getter)
        elif getter or setter:
            return self._make_decorator(getter, setter)
        else:
            return self._simulate
    
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
emptyblock._simulate = lambda: None
Block.emptyblock = emptyblock