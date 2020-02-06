# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:19:41 2019

@author: yoelr
"""
__all__ = ('design', 'add_design')


# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 21:18:50 2019

@author: yoelr
"""
from ..._stream import mol_flow_dim, mass_flow_dim, vol_flow_dim, _Q, DimensionError

__all__ = ('design',)

# %% Design Center class

def _design(self):
    D = self._Design
    for i, j in self._design_basis_: D[i] = j(self)

def _design_one(self):
    self._Design[self._design_name_] = self._design_func_()

class Basis:
    """Create a Basis object that defines a design basis. When called, it calls the factory function.
    
    **Parameters**
    
        **factory:** [function] Should return design basis function. Arguments should be 'units', 'N_ins', and '_N_outs', where:
            * **units:** [str] Units of measure
                
            * **N_ins:** [int] Number of input streams
            
            * **N_outs** [int] Number of output streams
        
        **cached:** [dict] Dictionary to cache basis functions
    
    """
    __slots__ = ('factory', 'cached')
    def __init__(self, factory):
        self.factory = factory
        self.cached = {}
        
    def __call__(self, units, N_ins, N_outs):
        """Return size/design function.
        
        **Parameters**
        
            **units:** [str] Units of measure
                
            **N_ins:** [int] Number of input streams
            
            **N_outs** [int] Number of output streams
        
        """
        args = (units, N_ins, N_outs)
        cached = self.cached
        if args in cached:
            func = cached[args]
        else:
            func = self.factory(units, N_ins, N_outs)
            cached[args] = func
        return func
    
    def __repr__(self):
        return f"<{type(self).__name__}: {self.factory.__name__}(units, N_ins, N_outs)>"    


class DesignCenter:
    """Create a DesignCenter object that manages all design basis. When called, it returns a Unit class decorator that adds a design item to the given Unit class."""
    
    def define(self, factory):
        """Define a new design basis.
        
        **Parameters**
        
            **factory:** [function] Should return design basis function. Arguments should be 'units', 'N_ins', and 'N_outs', where `units` are units of measure, and `N_ins` and `N_outs are the number of input and output streams respectively.
        
        .. Note::
            
            Design basis is registered with the name of the factory function.
        
        """
        name = factory.__name__
        if name in self:
            raise ValueError(f"basis '{name}' already implemented")
        self.__dict__[name] = basis = Basis(factory)
        return basis
    
    def __call__(self, name, units, fsize=None):    
        """Return a Unit class decorator that adds a size/design requirement to the class.
        
        **Parameters**
        
            **name:** Name of design item.
            
            **units:** Units of measure of design item.
            
            **fsize:** Should return design item given the Unit object. If None, defaults to function predefined for given name and units.
        
        """
        return lambda cls: self._add_design2cls(cls, name, units, fsize)
    
    def _add_design2cls(self, cls, name, units, fsize):
        """Add size/design requirement to class.
        
        **Parameters**
        
            **cls:** Unit class.
        
            **name:** Name of design item.
            
            **units:** Units of measure of design item.
            
            **fsize:** Should return design item given the Unit object. If None, defaults to function predefined for given name and units.
            
        **Examples**
        
            :doc:`Unit decorators`
        
        """
        f = fsize or design._get_fsize(name, units, cls._N_ins, cls._N_outs)
        
        # Make sure new _units dictionary is defined
        if not cls._units:
            cls._units = {}
        elif '_units' not in cls.__dict__:
            cls._units = cls._units.copy()
        
        # Make sure design basis is not defined
        if name in cls._units:
            raise RuntimeError(f"design basis '{name}' already defined in class")
        else:
            cls._units[name] = units
        
        # Add design basis
        if cls._design is _design:
            cls._design_basis_.append((name, f))
        elif cls._design is _design_one:
            cls._design_basis_ = [(cls._design_name_, cls._design_func_),
                                  (name, f)]
            cls._design = _design
            del cls._design_name_, cls._design_func_
        elif '_design' in cls.__dict__:
            raise RuntimeError("'_design' method already implemented")
        else:
            cls._design_name_ = name
            cls._design_func_ = f
            cls._design = _design_one
        
        return cls
    
    def _get_fsize(self, name, units, N_ins, N_outs):
        """Return size/design function.
        
        **Parameters**
        
            **name:** [str] name for design basis
            
            **units:** [str] Units of measure
                
            **N_ins:** [int] Number of input streams
            
            **N_outs** [int] Number of output streams
        
        """
        try:
            basis = getattr(self, name)
        except:
            raise ValueError(f"unknown basis '{name}', basis must be one of the following: {', '.join(self)}")
        return basis(units, N_ins, N_outs)

    def __getattr__(self, name):
        return object.__getattribute__(self, name.replace(' ', '_').casefold())

    def __setattr__(self, name, basis):
        raise AttributeError(f"can't set attribute")

    def __contains__(self, basis):
        return basis in self.__dict__
    
    def __iter__(self):
        yield from self.__dict__
    
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join(self)}>"

# %% Design factories
  
design = DesignCenter() #: Used to decorate classes with new design item

@design.define
def flow_rate(units, N_ins, N_outs):
    q = _Q(1, units)
    dim = q.dimensionality
    if dim == mol_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].molnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].molnet
        else: 
            func = lambda self: self._molnet_in
        ubase = 'kmol/hr'
    elif dim == mass_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].massnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].massnet
        else: 
            func = lambda self: self._massnet_in
        ubase = 'kg/hr'
    elif dim == vol_flow_dim:
        if N_ins == 1:
            func = lambda self: self._ins[0].volnet
        elif N_outs == 1:
            func = lambda self: self._outs[0].volnet
        else: 
            func = lambda self: self._volnet_in
        ubase = 'm3/hr'
    else:
        raise DimensionError(f"dimensions for flow units must be in molar, mass or volumetric flow rates, not '{dim}'")
    factor = 1/q.to(ubase).magnitude
    if factor == 1:
        return func
    else:
        return lambda self: factor*func(self)

@design.define
def duty(units, N_ins, N_outs):
    if not (N_ins == N_outs == 1):
        raise ValueError(f"number of input and output streams must be 1 for selected basis")
    factor = _Q(1, 'kJ/hr').to(units).magnitude
    if factor == 1:
        return lambda self: self._outs[0].H - self._ins[0].H
    else:
        def find_duty(self):
            self._duty_kJ_mol_ = duty = self._outs[0].H - self._ins[0].H
            return factor*duty
        return find_duty

@design.define
def dry_flow_rate(units, N_ins, N_outs):
    if not (N_ins == 1):
        raise ValueError(f"number of input streams must be 1 for selected basis")
    if units != 'kg/hr':
        raise ValueError(f"units must be in kg/hr for selected basis")
    def dry_flow(self):
        feed = self._ins[0]
        return feed.massnet - feed.mass[feed.index('7732-18-5')]
    return dry_flow

del flow_rate, duty, dry_flow_rate



