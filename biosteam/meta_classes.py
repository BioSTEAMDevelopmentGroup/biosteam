# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 13:13:41 2018

@author: Guest Group
"""
from .exceptions import metaError
from .unit import metaUnit

# %% For making Units final


class metaFinal(metaUnit):
    """Make Unit class final/unable to be inherited."""
    def __new__(mcl, name, bases, dct):
        # Make this class a final class
        for b in bases:
            if isinstance(b, mcl):
                raise metaError(
                    f"Cannot inherit from {b}. Instances of {mcl.__name__} cannot be inherited.")

        return super().__new__(mcl, name, bases, dct)


# %% Wrapping metaclass

def wrap_wrapper(func, new_func, order):
    if order == 'before':
        f1, f2 = new_func, func
    elif order == 'after':
        f2, f1 = new_func, func
    else:
        raise metaError(
            f"Invalid order. Only either 'before' or 'after' are valid, not {order}")

    def wrapper(self):
        f1(self)
        f2(self)

    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


class metaWrap(metaFinal):
    """Unit metaclass that wraps declared methods with superclass methods. Also adds new key word arguments specified in kwargs.

    **Parameters** (Class definitions)

         **order:** [str] if 'before', the declared methods will run before the superclass method definitions. If 'after', the run declared methods run after the superclass method definitions.

        **Optional**

             **kwargs:** [dict] New keyword arguments
    
             **setup():** Create components and cached data
    
             **run():** Run rigorous simulation
    
             **operation():** Find operation requirements
    
             **design():** Find design requirements
    
             **cost():** Find capital and annual cost

    .. Note::

         * Classes from metaWrap are end-of-the-line/final classes. No further inheritance can be made. 
    """
    def __new__(mcl, name, bases, dct):
        # Make sure an order was given
        order = dct.get('order')
        if not order:
            raise metaError("Missing class definition 'order'")

        # Make sure only one base is given
        if len(bases) != 1:
            raise metaError("One and only one superclass must be given")
        base = bases[0]

        # Include additional key word arguments
        dct['kwargs'].update(base.kwargs)

        # Wrap old methods with new methods
        for key in ('setup', 'run', 'operation', 'design', 'cost'):
            new_method = dct.get(key)
            if new_method:
                dct[key] = wrap_wrapper(getattr(base, key), new_method, order)

        # Return metaUnit class
        return super().__new__(name, bases, dct)
