# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
import pandas as pd
from warnings import warn
from ._graphics import UnitGraphics, box_graphics
from ._heat_utility import HeatUtility
from .utils import AbstractMethod, format_title, static, piping
from ._power_utility import PowerUtility
from .exceptions import UnitInheritanceError
from thermosteam.units_of_measure import convert
from copy import copy
import biosteam as bst
from thermosteam import Stream
from .utils.piping import InletPort, OutletPort, MissingStream, ignore_docking_warnings
from ._network import repr_ins_and_outs
from thermosteam.utils.decorators import registered, thermo_user
from thermosteam.utils import extended_signature
from thermosteam.base import display_asfunctor
from flexsolve import IQ_interpolation
from numpy.typing import NDArray
from typing import Callable, Optional, TYPE_CHECKING, Sequence, Iterable, Tuple, Any, Container
from thermosteam.base.sparse import SparseVector, sum_sparse_vectors
import numpy as np
import thermosteam as tmo
from thermosteam import EquationNode, VariableNode
if TYPE_CHECKING: 
    System = bst._recycle_system
    HXutility = bst.HXutility
    UtilityAgent = bst.UtilityAgent

streams = Optional[Sequence[Stream|str]]

__all__ = ('Unit', 'design_input', 'update_design_input', 'ProcessSpecification', 'BoundedNumericalSpecification')

int_types = (int, np.int32)

# %% Typing

# from typing import Collection, Union, Annotated
# streams = Union[Collection[Union[Stream, str, None]], Union[Stream, str, None]]
# stream = Union[Annotated[Union[Stream, str, None], 1], Union[Stream, str, None]]
# stream_sequence = Collection[Union[Stream, str, None]]

# %% Path utilities

def add_path_segment(start, end, path, ignored):
    fill_path_segment(start, path, end, set(), ignored)

def fill_path_segment(start, path, end, previous_units, ignored):
    if start is end: return path
    if start not in previous_units: 
        if start not in ignored: path.append(start)
        previous_units.add(start)
        success = False
        for outlet in start._outs:
            start = outlet._sink
            if not start: continue
            path_segment = fill_path_segment(start, [], end, previous_units, ignored)
            if path_segment is not None: 
                path.extend(path_segment)
                success = True
        if success: return path

# %% Auxiliary piping

def superposition_property(name):
    @property
    def p(self):
        return getattr(self.port.get_stream(), name)
    @p.setter
    def p(self, value):
        setattr(self.port.get_stream(), name, value)
        
    return p

def _superposition(cls, parent, port):
    # Stream slots for connectivity that are shared by default with superposition streams.
    _superposition_stream_slots = ('_ID', '_source', '_sink', '_thermo', 'port')
    
    # Stream methods and attributes that act on the superposition stream itself 
    # (not shared with the actual stream).
    _nonsuperposition_stream_names = (
        'line', 'feed_priorities', 'get_feed_priority', 'set_feed_priority', 
        'source', 'sink', 'isfeed', 'isproduct', 'disconnect_source', 
        'disconnect_sink', 'disconnect', '__sub__', '__rsub__', '__pow__', 
        '__rpow__', 'get_connection', '_basic_info', '_source_info',
        'registry', 'ticket_name', 'ticket_numbers', 'autonumber', 
        'unregistered_ticket_number', '_take_ticket', '_register', 
        'register_alias', 'ID', '__repr__', '__str__',
    )
    
    excluded = set([*cls.__dict__, *_nonsuperposition_stream_names, port, '_' + port, 'port', '_F_node', '_E_node', '_T_node'])
    for name in (*parent.__dict__, *_superposition_stream_slots):
        if name in excluded: continue
        setattr(cls, name, superposition_property(name))
    return cls

def superposition(parent, port):
    return lambda cls: _superposition(cls, parent, port)

@superposition(Stream, 'sink')
class SuperpositionInlet(Stream):
    """Stream that shares data with an inlet of a unit operation, but can be 
    connected to auxiliary units without being disconnected from the unit."""
    __slots__ = ()
    def __init__(self, port, sink=None):
        self.port = port
        self._sink = sink


@superposition(Stream, 'source')
class SuperpositionOutlet(Stream):
    """Stream that shares data with an outlet of a unit operation, but can be 
    connected to auxiliary units without being disconnected from the unit."""
    __slots__ = ()
    def __init__(self, port, source=None):
        self.port = port
        self._source = source

# %% Defining design inputs

def design_input(
        name: Optional[str]=None,
        units: Optional[str]=None,
        *,
        default: Optional[float|str]=None,
        lb: Optional[float]=None,
        ub: Optional[float]=None,
        catalog: Optional[Container[str]]=None,
        purchase_item: Optional[str]=None,
        catalog_type: Optional[str]=None,
        fget: Optional[Callable]=None,
        fset: Optional[Callable]=None,
        fdel: Optional[Callable]=None,
        doc: Optional[str]=None,
    ):
    """
    Return a design input descriptor. Design inputs behave like
    property objects (and are documented as such by Sphinx), but also
    validate values and are reported in :meth:`Unit.results`.

    Parameters
    ----------
    name :
        Name of design input. Defaults to the attribute name.
    units :
        Units of measure for numerical design input.
    default :
        Default value applied when design input is set to None.
    lb :
        Lower bound for numerical design input.
    ub :
        Upper bound for numerical design input.
    catalog :
        Valid options for categorical design input. If a dictionary,
        values are the cost factors of each option.
    purchase_item :
        Name of purchase cost item to apply cost factors.
    catalog_type :
        Type of cost factor: 'Material', 'Design', or 'Module'.
    fget :
        Getter. Defaults to returning the value stored.
    fset :
        Setter (called after validation). Defaults to storing the value.
    fdel :
        Deleter.
    doc :
        Docstring. Defaults to getter's docstring or "<name> [<units>].", if none available.

    Examples
    --------
    >>> import biosteam as bst
    >>> class Vessel(bst.Unit):
    ...     tau = bst.design_input('Residence time', 'hr', lb=0)
    ...     orientation = bst.design_input(catalog={'Vertical', 'Horizontal'}, default='Vertical')
    ...     open_roof = bst.design_input(default=False)
    ...     N_vessels: int = bst.design_input(default=1, lb=0)
    ...     
    ...     @bst.design_input(units='m3', lb=0)
    ...     def volume(self):
    ...         return self.tau * self.ins[0].F_vol
    ...     @volume.setter
    ...     def volume(self, volume):
    ...         self.tau = volume / self.ins[0].F_vol
    ...

    >>> Vessel.tau
    NumericalDesignInput('Residence time', 'hr', lb=0)
    >>> Vessel.orientation
    CategoricalDesignInput('Orientation', default='Vertical', catalog={'Horizontal', 'Vertical'})
    >>> Vessel.volume
    NumericalDesignInput('Volume', 'm3', lb=0)
    >>> Vessel.open_roof
    LogicalDesignInput('Open roof', default=False)
    >>> Vessel.N_vessels
    NumericalDesignInput('N vessels', default=1)
    
    """
    if catalog is None:
        if catalog_type is not None:
            raise ValueError('cannot define catalog_type without catalog')
        if purchase_item is not None:
            raise ValueError('cannot define purchase_item without catalog')
        if isinstance(default, bool):
            if units is not None:
                raise ValueError('logical design input cannot have units of measure')
            if lb is not None:
                raise ValueError('logical design input cannot have a lower bound (lb)')
            if ub is not None:
                raise ValueError('logical design input cannot have an upper bound (ub)')
            return LogicalDesignInput(name, default, fget, fset, fdel, doc)
        else:
            if units is None: units = ''
            return NumericalDesignInput(name, units, default, lb, ub, fget, fset, fdel, doc)
    else:
        N_none = (catalog_type is None) + (purchase_item is None)
        if N_none == 0:
            catalog_type = catalog_type.capitalize()
            if catalog_type == 'Material':
                catalog_type = 'F_M'
            elif catalog_type == 'Design':
                catalog_type = 'F_D'
            elif catalog_type == 'Module':
                catalog_type = 'F_BM'
            else:
                raise ValueError(f"{catalog_type!r} is not a recognised catalog_type ; only 'Material', 'Design', or 'Module' are valid")
        elif N_none == 2:
            purchase_item = None
            catalog_type = None
        else:
            raise ValueError('must define both catalog_type and purchase_item, or none of them; not either')
        if units is not None:
            raise ValueError('units cannot be defined together with catalog')
        if lb is not None:
            raise ValueError('lower bound (lb) cannot be defined together with catalog')
        if ub is not None:
            raise ValueError('upper bound (ub) cannot be defined together with catalog')
        return CategoricalDesignInput(name, default, catalog, purchase_item, catalog_type, fget, fset, fdel, doc)

class update_design_input:
    __slots__ = ('updates',)
    def __init__(self, **updates):
        self.updates = updates
    
    def __set_name__(self, owner, name):
        inherited = getattr(super(owner, owner), name, None)
        if not isinstance(inherited, DesignInput):
            raise RuntimeError(
                f'{owner.__name__}.{name} is not a design_input'
            )
        new = inherited.copy(**self.updates)
        new.__set_name__(owner, name)
        setattr(owner, name, new)
    
    def __repr__(self):
        kwargs = ', '.join([f'{i}={j!r}' for i, j in self.updates.items()])
        return f"{type(self).__name__}({kwargs})"

def _get_class_annotation(owner, key):
    try:
        import annotationlib
    except ImportError: # Python < 3.14
        annotations = owner.__dict__.get('__annotations__', {})
    else:
        try:
            annotations = annotationlib.get_annotations(
                owner, format=annotationlib.Format.FORWARDREF
            )
        except Exception:
            return None
    return annotations.get(key)


class DesignInput:
    
    def copy(self, **kwargs):
        for i in ('name', 'default', 'fget', 'fset', 'fdel', 'doc', *self._new_attributes):
            if i not in kwargs:
                kwargs[i] = getattr(self, i)
        copy = type(self)(**kwargs)
        if hasattr(self, '_return_type'):
            copy._default_return_type = self._return_type
        if hasattr(self, '__name__'):
            copy.__name__ = self.__name__
        return copy

    def __set_name__(self, owner, name):
        self.__name__ = name
        if self.name is None: self.name = name.replace('_', ' ').capitalize()
        fset = self.fset
        fget = self.fget
        N_none = (fget is None) + (fset is None)
        
        # Set return type for validation
        if N_none == 1:
            raise RuntimeError('design inputs must implement both getter and setter methods, not just one')
        elif N_none == 0:
            annotations = fget.__annotations__
            if 'return' in annotations:
                self._return_type = annotations['return']
            else:
                annotations['return'] = self._return_type = self._default_return_type
        else:
            self._return_type = _get_class_annotation(owner, name) or self._default_return_type
            
        # Validate default
        default = self.default
        if default is not None: 
            self.default = default = self.validate(default)
        
        # Setup docstring
        doc = self.doc
        if doc is None:
            if fget is None or (doc:=fget.__doc__) is None:
                units = self.units
                if default is None:
                    default_note = ''
                elif units:
                    default_note = f' Defaults to {default!r} [{units}].'
                else:
                    default_note = f' Defaults to {default!r}.'
                if units:
                    doc = f"{self.name} [{units}].{default_note}" 
                else:
                    doc = f"{self.name}.{default_note}"
            self.__doc__ = doc
            
        # Default getter/setter
        if N_none == 2: 
            # Document like a property with a getter named after the attribute
            attr = '_' + name
            self.fget = eval(f'lambda self: self.{attr}')
            # fget.__name__ = name
            # fget.__qualname__ = f"{owner.__qualname__}.{name}"
            # fget.__module__ = owner.__module__
            # fget.__doc__ = doc
            # fget.__annotations__ = {'return': self._return_type}
            self.fset = eval(f'lambda unit, value: setattr(unit, {attr!r}, value)')

    def __call__(self, fget):
        # Used as a decorator right after creation (e.g., @design_input(...)),
        # so the getter is set in place instead of creating a copy
        self.fget = fget
        return self

    def getter(self, fget):
        if hasattr(self, '__name__'): return self.copy(fget=fget)
        self.fget = fget
        return self

    def setter(self, fset):
        if hasattr(self, '__name__'): return self.copy(fset=fset)
        self.fset = fset
        return self

    def deleter(self, fdel):
        if hasattr(self, '__name__'): return self.copy(fdel=fdel)
        self.fdel = fdel
        return self

    def __get__(self, instance, owner):
        if instance is None: return getattr(owner, self.__name__)
        return self.fget(instance)

    def __delete__(self, unit):
        if self.fdel is None: raise AttributeError(f'cannot delete design input {self.__name__!r}')
        self.fdel(unit)


class NumericalDesignInput(DesignInput):
    _new_attributes = ('units', 'lb', 'ub')
    _default_return_type = float

    def __init__(self, name, units, default, lb, ub, fget, fset, fdel, doc):
        self.name = name
        self.units = units
        self.default = default
        self.lb = lb
        self.ub = ub
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.doc = doc
    
    def validate(self, value):
        new_value = self._return_type(value)
        if new_value == value:
            value = new_value
            if self.lb is not None and value < self.lb:
                units = self.units
                units = ' ' + units if units else ''
                try:
                    msg = f"{self.__name__} ({value:.4g}{units}) cannot be below {self.lb:.4g}{units}"
                except:  # Handle format errors
                    msg = f"{self.__name__} ({value}{units}) cannot be below {self.lb}{units}"
            elif self.ub is not None and value > self.ub:
                units = self.units
                units = ' ' + units if units else ''
                try:
                    msg = f"{self.__name__} ({value:.4g}{units}) cannot be above {self.ub:.4g}{units}"
                except:  # Handle format errors
                    msg = f"{self.__name__} ({value}{units}) cannot be above {self.ub}{units}"
            else:
                return value
        else:
            msg = f'{self.__name__} must be a {self._return_type.__name__}; {type(value).__name__} object is invalid'
        raise AttributeError(msg)

    def __set__(self, unit, value):
        if value is not None: self.validate(value)
        self.fset(unit, value)

    def __repr__(self):
        default = '' if self.default is None else f', default={self.default!r}'
        units = f', {self.units!r}' if self.units else ''
        lb = '' if self.lb is None else f', lb={self.lb}'
        ub = '' if self.ub is None else f', ub={self.ub}'
        return f"{type(self).__name__}({self.name!r}{units}{default}{lb}{ub})"


class CategoricalDesignInput(DesignInput):
    _new_attributes = ('catalog', 'purchase_item', 'catalog_type')
    units = ''
    _default_return_type = str

    def __init__(self, name, default, catalog, purchase_item, catalog_type, fget, fset, fdel, doc):
        self.name = name
        self.default = default
        self.catalog = catalog
        self.purchase_item = purchase_item
        self.catalog_type = catalog_type
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.doc = doc

    def validate(self, value):
        catalog = self.catalog
        if value not in catalog:
            try: 
                # Assume value is a string
                value = value.capitalize()
            except: 
                # If its not (e.g., could be an int), that is fine
                pass
            else: 
                if value in catalog: return value
            N = len(catalog)
            if N > 2:
                *names, last = list(catalog)
                names = ', '.join([repr(i) for i in names]) + ', and ' + repr(last)
                msg = f"{value!r} is not a recognised {self.__name__}; only {names} are valid"
            elif N == 2:
                left, right = catalog
                msg = f"{value!r} is not a recognised {self.__name__}; only {left!r} and {right!r} are valid"
            elif N == 1:
                name = list(catalog)[0]
                msg = f"{value!r} is not a recognised {self.__name__}; only {name} is valid"
            else:
                msg = f"no valid {self.__name__} available"
        else:
            return value
        raise AttributeError(msg)

    def __set__(self, unit, value):
        if value is None:
            if self.purchase_item is not None:
                getattr(unit, self.catalog_type).pop(self.purchase_item, None)
        else:
            value = self.validate(value)
            if self.purchase_item is not None:
                getattr(unit, self.catalog_type)[self.purchase_item] = self.catalog[value]
        self.fset(unit, value)

    def __repr__(self):
        catalog = self.catalog
        if isinstance(catalog, (set, frozenset)):
            catalog = ', catalog={' + str(sorted(catalog))[1:-1] + '}'
        elif catalog is None:
            catalog = ''
        else:
            catalog = f', catalog={catalog}'
        catalog = '' if self.catalog is None else f', catalog={self.catalog}'
        purchase_item = '' if self.purchase_item is None else f', purchase_item={self.purchase_item}'
        catalog_type = '' if self.catalog_type is None else f', catalog_type={self.catalog_type}'
        default = '' if self.default is None else f', default={self.default!r}'
        return f"{type(self).__name__}({self.name!r}{default}{catalog}{purchase_item}{catalog_type})"


class LogicalDesignInput(DesignInput):
    _new_attributes = ()
    _default_return_type = bool
    units = ''

    def __init__(self, name, default, fget, fset, fdel, doc):
        self.name = name
        self.default = default
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.doc = doc

    def validate(self, value):
        if isinstance(value, bool):
            return value
        else:
            raise AttributeError(f"{self.__name__} must be either True or False; {type(value).__name__} object is invalid")

    def __set__(self, unit, value):
        value = self.validate(value)
        self.fset(unit, value)

    def __repr__(self):
        default = '' if self.default is None else f', default={self.default!r}'
        return f"{type(self).__name__}({self.name!r}{default})"
    

# %% Unit Operation

@thermo_user
@registered(ticket_name='U')
class Unit:
    """
    Abstract class for Unit objects. Subclasses must define
    :attr:`~Unit._run`, :attr:`~Unit._design` and :attr:`~Unit._cost` methods to 
    estimate stream outlets of a Unit and find design and cost information.  
    
    For advanced phenomena based simulation, subclasses may define 
    :attr:`~Unit._update_nonlinearities`,
    :attr:`~Unit._create_bulk_balance_equations`,
    :attr:`~Unit._create_material_balance_equations`,
    :attr:`~Unit._create_energy_balance_equations`,
    :attr:`~Unit._update_energy_coefficient`,
    and :attr:`~Unit._update_variable` methods.
    
    Parameters
    ----------
    ID :
        A unique identification. If ID is None, the unit will not be
        registered in flowsheet. The ID defaults to the variable assigned in 
        the source code. If no variable is assigned, a unique ID will be 
        chosen.
    ins :
        Inlet streams or IDs to initialize inlet streams.
        If empty tuple, streams with default IDs will be created.
        By default, streams will be missing.
    outs : 
        Outlet streams or IDs to initialize outlet streams.
        By default, streams with unique IDs will be created.
        If None, streams will be missing.
    thermo : 
        Thermo object to initialize inlet and outlet streams. Defaults to
        :meth:`settings.thermo <thermosteam._settings.ProcessSettings.thermo>`.
    
    Notes
    -----
    The free on board (f.o.b.) purchase costs and installed equipment costs 
    (i.e. bare-module cost) for each item in the :attr:`~Unit.baseline_purchase_costs` 
    dictionary and in auxiliary units are automatically added to the 
    :attr:`~Unit.purchase_costs` and :attr:`~Unit.installed_costs` dictionaries. 
    
    As explained in [1]_, the f.o.b. purchase cost is given by:
    
    .. math::
       
       C_{P} = C_{Pb}F_{D}F_{P}F_{M}
    
    And the installed equipment cost is given by:
    
    .. math::
       
       C_{BM} = C_{Pb} (F_{BM} + F_{D}F_{P}F_{M} - 1)
    
    Where:
        * :math:`C_{Pb}`: Baseline purchase cost.
        * :math:`F_{BM}`: Bare module factor.
        * :math:`F_{D}`: Design factor.
        * :math:`F_{P}`: Pressure factor.
        * :math:`F_{M}`: Material factor.
    
    Values for the bare-module, design, pressure, and material factors of 
    each equipment should be stored in the :attr:`~Unit.F_BM`, :attr:`~Unit.F_D`, 
    :attr:`~Unit.F_P`, and :attr:`~Unit.F_M` dictionaries.
    
    Examples
    --------
    :doc:`../tutorial/Creating_a_Unit`
    
    :doc:`../tutorial/-pipe-_notation`
    
    :doc:`../tutorial/Inheriting_from_Unit`
    
    References
    ----------
    .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley. Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    
    """ 
    max_parallel_units = int(10e3)
    _design_inputs = []
    
    ### Abstract Attributes ###
    
    #: **class-attribute** Units of measure for :attr:`~Unit.design_results` dictionary.
    _units: dict[str, str] = {}
    
    #: **class-attribute** Expected number of inlet streams. Defaults to 1.
    _N_ins: int = 1  
    
    #: **class-attribute** Expected number of outlet streams. Defaults to 1
    _N_outs: int = 1
    
    #: **class-attribute** Whether the number of streams in :attr:`~Unit.ins` is fixed.
    _ins_size_is_fixed: bool = True
    
    #: **class-attribute** Whether the number of streams in :attr:`~Unit.outs` is fixed.
    _outs_size_is_fixed: bool = True
    
    #: **class-attribute** Name denoting the type of Unit class. Defaults to the class
    #: name of the first child class
    line: str = 'Unit'
    
    #: **class-attribute** Default bare-module factors for each purchase cost item.
    #: Items in this dictionary are copied to the :attr:`~Unit.F_BM` attribute during 
    #: initialization.
    _F_BM_default: dict[str, float] = {}
    
    #: **class-attribute** Lifetime of equipment. Defaults to lifetime of
    #: production venture. Use an integer to specify the lifetime for all
    #: items in the unit purchase costs. Use a dictionary to specify the 
    #: lifetime of each purchase cost item.
    _default_equipment_lifetime: int|dict[str, int] = {}
    
    #: **class-attribute** Name of attributes that are auxiliary units. These units
    #: will be accounted for in the purchase and installed equipment costs
    #: without having to add these costs in the :attr:`~Unit.baseline_purchase_costs` dictionary.
    #: Heat and power utilities are also automatically accounted for.
    auxiliary_unit_names: tuple[str, ...] = ()
    
    #: **class-attribute** Cost items that need to be summed across operation modes for 
    #: flexible operation (e.g., filtration membranes).
    _materials_and_maintenance: frozenset[str] = frozenset()
    
    #: **class-attribute** Settings for diagram representation. Defaults to a 
    #: box with the same number of inlet and outlet edges as :attr:`~Unit._N_ins` 
    #: and :attr:`~Unit._N_outs`.
    _graphics: UnitGraphics = box_graphics
    
    #: **class-attribute** Whether to skip detailed simulation when inlet 
    #: streams are empty. If inlets are empty and this flag is True,
    #: detailed mass and energy balance, design, and costing algorithms are skipped
    #: and all outlet streams are emptied.
    _skip_simulation_when_inlets_are_empty = False
    
    #: **class-attribute** Whether to link inlet and outlet streams.
    _link_streams: bool = False

    #: Whether unit leverages information from other units within the system 
    #: regardless of whether or not they are directly connected.
    _universal = False
    
    #: Whether streams do not mix, but do interact.
    _interaction = False
    
    #: **class-attribute** Index for auxiliary inlets to parent unit for graphviz diagram settings.
    _auxin_index = {}

    #: **class-attribute** Index for auxiliary outlets to parent unit for graphviz diagram settings.
    _auxout_index = {}
    
    #: **class-attribute** Used for piping warnings.
    _stacklevel: int = 6

    ### Abstract methods ###
    
    #: Create auxiliary components.
    _load_components = AbstractMethod
    
    #: Add design requirements to the :attr:`~Unit.design_results` dictionary.
    _design = AbstractMethod
    
    #: Add itemized purchase costs to the :attr:`~Unit.baseline_purchase_costs` dictionary.
    _cost = AbstractMethod

    #: Add embodied emissions (e.g., unit construction) in LCA.
    _lca = AbstractMethod
    
    #: Initialize unit operation with key-word arguments.
    _init = AbstractMethod
    
    #: Run mass and energy balances and update outlet streams (without user-defined specifications).
    _run = AbstractMethod
    
    def __init_subclass__(cls,
            isabstract=False,
            new_graphics=True,
            does_nothing=None,
        ):
        dct = cls.__dict__
        # Collect design inputs through the class history (including mixins);
        # attributes that override an inherited design input remove it.
        design_inputs = {}
        for base in reversed(cls.__mro__):
            for key, obj in base.__dict__.items():
                if isinstance(obj, DesignInput):
                    design_inputs[key] = obj
                elif key in design_inputs:
                    del design_inputs[key]
        design_inputs = [i for i in design_inputs.values()]
        cls._design_inputs = [i.__name__ for i in design_inputs] # Preserve order
        
        # Extend annotation
        if '__init__' in dct and '_init' not in dct:
            init = dct['__init__']
            if hasattr(init, 'extension'): 
                cls._init = init.extension
        elif dct.get('_init'):
            _init = dct['_init']
            cls.__init__ = extended_signature(cls.__init__, _init, design_inputs)
        else:
            for i in dct.values(): 
                if not isinstance(i, DesignInput): continue
                cls.__init__ = extended_signature(cls.__init__, design_inputs=design_inputs)
                break
                
        # Standard annotation and adjustment to stacklevel
        if '__init__' in dct or '_init' in dct:
            init = dct['__init__']
            annotations = init.__annotations__
            for i in ('ins', 'outs'):
                if i not in annotations: annotations[i] = streams
            if '_stacklevel' not in dct: cls._stacklevel += 1
        
        if does_nothing: return 
        dct = cls.__dict__
        if 'run' in dct:
            raise UnitInheritanceError(
                 "the 'run' method cannot be overridden; implement `_run` instead"
            )
        if '_N_heat_utilities' in dct:
            warn("'_N_heat_utilities' class attribute is scheduled for deprecation; "
                 "use the `add_heat_utility` method instead",
                 DeprecationWarning, stacklevel=2)
        if 'line' not in dct:
            cls.line = format_title(cls.__name__)
        if 'ticket_name' not in dct:
            line = cls.line.lower()
            if 'centrifuge' in line: cls.ticket_name = 'C'
            elif 'distillation' in line: cls.ticket_name = 'D'
            elif 'evaporator' in line: cls.ticket_name = 'E'
            elif 'flash' in line: cls.ticket_name = 'F'
            elif ('cooler' in line 
                  or 'condenser' in line
                  or 'heater' in line 
                  or 'boiler' in line
                  or 'heat exchanger' in line): cls.ticket_name = 'H'
            elif 'compressor' in line: cls.ticket_name = 'K'
            elif 'turbine' in line: cls.ticket_name = 'Ʞ'
            elif 'mixer' in line: cls.ticket_name = 'M'
            elif 'pump' in line: cls.ticket_name = 'P'
            elif 'reactor' in line or 'digestion' in line or 'ferment' in line: cls.ticket_name = 'R'
            elif 'split' in line: cls.ticket_name = 'S'
            elif 'tank' in line: cls.ticket_name = 'T'
            elif 'junction' == line: cls.ticket_name = 'J'
            elif 'specification' in line: cls.ticket_name = 'PS'
            elif 'valve' in line: cls.ticket_name = 'V'
            else: cls.ticket_name = 'U'
        if '_graphics' not in dct and new_graphics:
            # Set new graphics for specified line
            cls._graphics = UnitGraphics.box(cls._N_ins, cls._N_outs)
        if not isabstract:
            if cls._F_BM_default is Unit._F_BM_default:
                cls._F_BM_default = {}
            if cls._default_equipment_lifetime is Unit._default_equipment_lifetime: 
                cls._default_equipment_lifetime = {}
            if cls._units is Unit._units: 
                cls._units = {}
            if not cls._run:
                if cls._N_ins == 1 and cls._N_outs == 1:
                    static(cls)
                elif (cls._N_ins or cls._N_outs):
                    raise UnitInheritanceError(
                        "Unit subclass with multiple inlet or outlet streams "
                        "must implement a '_run' method unless the "
                        "'isabstract' keyword argument is True"
                    )
        name = cls.__name__
        if hasattr(bst, 'units') and hasattr(bst, 'wastewater') and hasattr(bst, 'facilities'):
            # Add 3rd party unit to biosteam module for convenience
            if isinstance(cls, bst.Facility):
                if name not in bst.facilities.__dict__:
                    bst.facilities.__dict__[name] = cls
            elif name not in bst.units.__dict__:
                bst.units.__dict__[name] = cls
            if name not in bst.__dict__:
                bst.__dict__[name] = cls
    
    def _collect_variables(self, key):
        self._system_collect_variables((key, self.ID))
    
    def get_E_node(self, stream):
        return getattr(self, 'E_node', None)
    
    def initialize_equation_nodes(self):
        for eq in self.equation_node_names: getattr(self, f'initialize_{eq}')()
        
    @property
    def variable_nodes(self):
        return set(sum([i.variable_nodes for i in self.equation_nodes], []))
    
    def create_specification_nodes(self, stream_ref):
        self.material_balance_specifications_nodes = specification_nodes = []
        self.specification_streams = specification_streams = []
        self.specification_equations = specification_equations = []
        for s in (*self.ins, *self.outs):
            for i, f in enumerate(s.material_equations):
                node = EquationNode(f"{self.node_tag}.material_balance_specifications_nodes[{i}]")
                specification_nodes.append(node) 
                specification_equations.append(f)
                dct, _ = f()
                streams = [stream_ref[i.material_reference] for i in dct]
                inputs = []
                outputs = []
                for s in streams:
                    isfeed = s.isfeed() or s.source._recycle_system is not self._recycle_system
                    if isfeed:   
                        sink = s.sink
                        if not hasattr(s, '_F_node'):
                            s._F_node = VariableNode(
                                f"{sink.node_tag}.ins[{sink.ins.index(s)}].F",
                                lambda: s.mol.to_array(),
                            )
                        outputs.append(s.F_node)
                    else:
                        inputs.append(s.F_node)
                specification_streams.extend(streams)
                node.set_equations(inputs=inputs, outputs=outputs)
    
    def _collect_specification_errors(self):
        results = []
        for node, f in zip(self.material_balance_specifications_nodes, 
                           self.specification_equations):
            dct, rhs = f()
            lhs = sum([s.mol * coef for s, coef in dct.items()])
            error = np.abs(rhs - lhs).sum()
            results.append((node.name, error))
        return results
    
    def create_equation_nodes(self, stream_ref):
        eqs = self.equation_node_names
        equation_nodes = []
        for eq in eqs: 
            node = EquationNode(f"{self.node_tag}.{eq}")
            equation_nodes.append(node)
            setattr(self, eq, node)
        self.create_specification_nodes(stream_ref)
        self.equation_nodes = (*equation_nodes, *self.material_balance_specifications_nodes)
        if not hasattr(self, 'energy_balance_node'):
            self.energy_balance_node = None
    
    @property
    def node_tag(self):
        if self.owner is self:
            return self.ID
        else:
            return f"{self.owner.node_tag}.{self.ID}"
        
    def _update_nonlinearities(self):
        """
        Update phenomenological variables for phenomena-based simulation.
        """
        ins = self.ins
        outs = self.outs
        data = [i.get_data() for i in outs]
        Ts = [i.T for i in outs]
        Q = self.Hnet
        self._run()
        f = bst.PhasePartition.dmol_relaxation_factor
        old = self._dmol
        new = sum(
            [i.mol for i in outs],
            -sum([i.mol for i in ins], 0)
        )
        self._dmol = dmol = f * old + (1 - f) * new
        self._dmol[np.abs(dmol) < 1e-9] = 0.
        self._duty = self.Hnet
        Ts_new = [i.T for i in outs]
        if all([i == j for i, j in zip(Ts_new, Ts)]): # T constant
            self.specified_variable = 'T'
        elif self.Hnet == Q: # Q constant    
            self.specified_variable = 'Q'
        else:
            self.specified_variable = 'B'
        for i, j in zip(outs, data): i.set_data(j)
    
    def _simulation_error(self):
        ins = self.ins
        outs = self.outs
        new_ins = [i.copy() for i in ins]
        new_outs = [i.copy() for i in outs]
        streams = new_ins + new_outs
        Ts = np.array([i.T for i in streams])
        flows = np.array([i.mol for i in streams])
        self._ins = new_ins
        self._outs = new_outs
        self._run()
        Ts_new = np.array([i.T for i in streams])
        new_flows = np.array([i.mol for i in streams])
        self._ins = ins
        self._outs = outs
        return np.abs(new_flows - flows).sum(), np.abs(Ts_new - Ts).sum()
    
    def _update_energy_coefficient(self, stream, coefficients):
        """
        Update coefficient of a stream 
        for phenomena-based simulation and 
        return the reference value.
        """
        if self.specified_variable == 'Q': 
            C = stream.C
            coefficients[self, 'T'] = -C
            return -C * stream.T
        return 0
    
    def _update_variable(self, variable, value):
        """
        Update variable being solved in equations for 
        phenomena-based simulation.
        """
        if variable == 'T':
            for i in self.outs: i.T = value
        else:
            raise ValueError('invalid variable')
    
    def _create_material_balance_equations(self, composition_sensitive):
        """
        list[tuple[dict, array]] Create material balance equations for 
        phenomena-based simulation.
        """
        if self._link_streams: return []
        fresh_inlets, process_inlets, equations = self._begin_material_equations(composition_sensitive)
        outs = self.flat_outs
        N = self.chemicals.size
        ones = np.ones(N)
        predetermined_flow = SparseVector.from_dict(sum_sparse_vectors([i.mol for i in fresh_inlets]), size=N)
        try:
            dmol = self._dmol
            rhs = predetermined_flow + dmol
        except:
            rhs = predetermined_flow
        mol_total = sum([i.mol for i in outs])
        for s in outs:
            split = s.mol / mol_total
            minus_split = -split
            eq_outs = {}
            for i in process_inlets: eq_outs[i] = minus_split
            eq_outs[s] = ones
            equations.append(
                (eq_outs, split * rhs)
            )
        return equations
    
    def _create_energy_balance_equations(self):
        """
        list[tuple[dict, float]] Create energy departure equations for 
        phenomena-based simulation.
        """
        fresh_inlets, process_inlets, equations = self._begin_energy_equations()
        if self.specified_variable == 'Q':
            dmol = self._dmol
            if dmol.any(): 
                dH = self._duty + sum([i.Hf for i in self.outs]) - sum([i.Hf for i in self.ins]) 
            else:
                dH = self._duty
            coeff = {(self, 'T'): sum([i.C for i in self.outs])}
            for i, s in enumerate(self.outs): coeff[self, i] = s.h
            for i in process_inlets: 
                dH += i._update_energy_coefficient(coeff)
            for i in fresh_inlets:
                dH += i.H
            equations.append(
                (coeff, dH)
            )
            return equations
        else:
            return equations
    
    def _create_bulk_balance_equations(self):
        """
        list[tuple[dict, array]] Create bulk balance equations for 
        phenomena-based simulation.
        """
        if self._link_streams: return []
        fresh_inlets, process_inlets, equations = self._begin_bulk_equations()
        outs = self.flat_outs
        predetermined_flow = [i.F_mol for i in fresh_inlets]
        try:
            dF_mol = self._dmol.sum()
            rhs = predetermined_flow + dF_mol
        except:
            rhs = predetermined_flow
        F_mol_total = sum([i.F_mol for i in outs])
        for i, s in enumerate(outs):
            split = s.F_mol / F_mol_total
            minus_split = -split
            eq_outs = {}
            for i in process_inlets: eq_outs[i, 'F_mol'] = minus_split
            eq_outs[s, 'F_mol'] = 1
            equations.append(
                (eq_outs, split * rhs)
            )
        return equations
    
    def _begin_material_equations(self, composition_sensitive):
        inlets = self.ins
        fresh_inlets = []
        process_inlets = []
        material_equations = [f() for i in inlets for f in i.material_equations]
        isfeed = lambda s: s.isfeed() or s.source._recycle_system is not self._recycle_system
        if not material_equations:
            if composition_sensitive:
                for i in inlets: 
                    if (isfeed(i) or
                        not getattr(i.source, 'composition_sensitive', False)):
                        fresh_inlets.append(i)
                    elif len(i) > 1:
                        process_inlets.extend(i)
                    else:
                        process_inlets.append(i)                    
            else:
                for i in inlets: 
                    if isfeed(i):
                        fresh_inlets.append(i)
                    elif len(i) > 1:
                        process_inlets.extend(i)
                    else:
                        process_inlets.append(i)   
            equations = material_equations
        elif composition_sensitive:
            equations = []
            dependent_streams = []
            for eq in material_equations:
                dct, value = eq 
                for i in dct:
                    if (i.source and
                        i.source._recycle_system is self._recycle_system and 
                        not getattr(i.source, 'composition_sensitive', False)):
                        break
                else:
                    dependent_streams.extend([i.imol for i in dct])
                    equations.append(eq)
            dependent_streams = set(dependent_streams)
            for i in inlets: 
                if (i.imol not in dependent_streams and isfeed(i) or
                    not getattr(i.source, 'composition_sensitive', False)):
                    fresh_inlets.append(i)
                elif len(i) > 1:
                    process_inlets.extend(i)
                else:
                    process_inlets.append(i)   
        else:
            for i in inlets: 
                if (isfeed(i) and 
                    not i.material_equations):
                    fresh_inlets.append(i)
                elif len(i) > 1:
                    process_inlets.extend(i)
                else:
                    process_inlets.append(i)   
            equations = material_equations
        return fresh_inlets, process_inlets, equations
    
    def _begin_bulk_equations(self):
        inlets = self.ins
        fresh_inlets = []
        process_inlets = []
        equations = [f() for i in inlets for f in i.bulk_equations]
        isfeed = lambda s: s.isfeed() or s.source._recycle_system is not self._recycle_system
        if equations:
            for i in inlets: 
                if (isfeed(i) and 
                    not i.bulk_equations):
                    fresh_inlets.append(i)
                elif len(i) > 1:
                    process_inlets.extend(i)
                else:
                    process_inlets.append(i)  
        else:
             for i in inlets: 
                 if isfeed(i):
                     fresh_inlets.append(i)
                 elif len(i) > 1:
                     process_inlets.extend(i)
                 else:
                     process_inlets.append(i)
        return fresh_inlets, process_inlets, equations
    
    def _begin_energy_equations(self):
        inlets = self.ins
        fresh_inlets = []
        process_inlets = []
        equations = [f() for i in inlets for f in i.energy_equations]
        isfeed = lambda s: s.isfeed() or s.source._recycle_system is not self._recycle_system
        if equations:
            for i in inlets: 
                if (isfeed(i) and 
                    not i.energy_equations):
                    fresh_inlets.append(i)
                elif len(i) > 1:
                    process_inlets.extend(i)
                else:
                    process_inlets.append(i)  
        else:
             for i in inlets: 
                 if isfeed(i):
                     fresh_inlets.append(i)
                 elif len(i) > 1:
                     process_inlets.extend(i)
                 else:
                     process_inlets.append(i)
        return fresh_inlets, process_inlets, equations    
    
    def __init__(self,
            ID: Optional[str]='',
            ins: streams=None,
            outs: streams=(),
            thermo: Optional[tmo.Thermo]=None,
            **kwargs
        ):
        ### Initialize utilities
    
        #: All heat utilities associated to unit. Cooling and heating requirements 
        #: are stored here (including auxiliary requirements).
        self.heat_utilities: list[HeatUtility, ...] = [HeatUtility for i in range(getattr(self, '_N_heat_utilities', 0))]
        
        #: Electric utility associated to unit (including auxiliary requirements).
        self.power_utility: PowerUtility = PowerUtility()
    
        ### Initialize design/cost/LCA results
        
        try:
            #: All bare-module factors for each purchase cost. Defaults to values in 
            #: the class attribute :attr:`~Unit._F_BM_default`.
            self.F_BM: dict[str, float] = self._F_BM_default.copy()
        except AttributeError:
            self.F_BM = {}
        
        #: All design factors for each purchase cost item in :attr:`~Unit.baseline_purchase_costs`.
        self.F_D: dict[str, float] = {}
        
        #: All pressure factors for each purchase cost item in :attr:`~Unit.baseline_purchase_costs`.
        self.F_P: dict[str, float] = {}
        
        #: All material factors for each purchase cost item in :attr:`~Unit.baseline_purchase_costs`.
        self.F_M: dict[str, float] = {}
        
        #: All design requirements excluding utility requirements and detailed 
        #: auxiliary unit requirements.
        self.design_results: dict[str, object] = {}
        
        #: All baseline purchase costs without accounting for design, 
        #: pressure, and material factors.
        self.baseline_purchase_costs: dict[str, float] = {}
        
        #: Itemized purchase costs (including auxiliary units)
        #: accounting for design, pressure, and material factors (i.e., 
        #: :attr:`~Unit.F_D`, :attr:`~Unit.F_P`, :attr:`~Unit.F_M`).
        #: Items here are automatically updated at the end of unit simulation.
        self.purchase_costs: dict[str, float] = {}
        
        #: All installed costs accounting for bare module, design, 
        #: pressure, and material factors. Items here are automatically updated
        #: at the end of unit simulation.
        self.installed_costs: dict[str, float] = {}
        
        #: Indices of additional utilities given by inlet streams.
        self._inlet_utility_indices: dict[str, int] = {}
        
        #: Indices of additional utilities given by outlet streams.
        self._outlet_utility_indices: dict[str, int] = {}
        
        #: Indices of additional credits/fees given by inlet streams.
        self._inlet_cost_indices: dict[str, int] = {}
        
        #: Indices of additional credits/fees given by outlet streams.
        self._outlet_revenue_indices: dict[str, int] = {}
        
        #: Fractions of stream flow rates that receive fees or credits.
        self._flow_fractions: dict[str, int] = {}
        
        try:
            #: Lifetime of equipment. Defaults to values in the class attribute 
            #: :attr:`~Unit._default_equipment_lifetime`. Use an integer to specify the lifetime 
            #: for all items in the unit purchase costs. Use a dictionary to specify 
            #: the lifetime of each purchase cost item.
            self.equipment_lifetime: int|dict[str, int] = copy(self._default_equipment_lifetime)
        except AttributeError:
            self.equipment_lifetime = {}
        
        #: Name-number pairs of baseline purchase costs and auxiliary unit 
        #: operations in parallel. Use 'self' to refer to the main unit. Capital 
        #: and heat and power utilities in parallel will become proportional to this 
        #: value.
        self.parallel: dict[str, int] = {}
        
        #: Unit design decisions that must be solved to satisfy specifications.
        #: While adding responses is optional, simulations benefit from responses
        #: by being able to predict better guesses.
        self.responses: set[bst.GenericResponse] = set()
        
        self._utility_cost = None
        
        self._recycle_system = None
        
        if bst.settings.ID_magic and ID == '': ID = bst.utils.infer_variable_assignment(self.__class__)
        
        self._isdynamic = False
        self._system = None
        self._register(ID)
        self._load_thermo(thermo)
    
        ### Initialize streams
        self.auxins = {} #: dict[int, stream] Auxiliary inlets by index.
        self.auxouts = {} #:  dict[int, stream] Auxiliary outlets by index.
        self._init_inlets(ins)
        self._init_outlets(outs)
    
        ### Initialize specification    
    
        #: All specification functions
        self._specifications: list[Callable] = []
        
        #: Whether to run mass and energy balance after calling
        #: specification functions
        self.run_after_specifications: bool = False 
        
        #: Whether to prioritize unit operation specification within recycle loop (if any).
        self.prioritize: bool = False
        
        #: Safety toggle to prevent infinite recursion.
        self._active_specifications: set[ProcessSpecification] = set()
        
        #: Auxiliary unit operation names.
        self.auxiliary_unit_names = list(self.auxiliary_unit_names)
        
        # Set defaults for design results.
        cls = type(self)
        for i in self._design_inputs:
            if i in kwargs: 
                setattr(self, i, kwargs.pop(i))
            else:
                setattr(self, i, getattr(cls, i).default)
        
        self._init(**kwargs)

        self._assert_compatible_property_package()
        
    __init__._init_inputs = frozenset()
        
    def _init_inlets(self, ins):
        self._ins = piping.Inlets(
            self, self._N_ins, ins, self._thermo, self._ins_size_is_fixed, self._stacklevel,
        )
        
    def _init_outlets(self, outs, **kwargs):
        self._outs = piping.Outlets(
            self, self._N_outs, outs, self._thermo, self._outs_size_is_fixed, self._stacklevel,
        )
        
    @property
    def ins(self) -> Sequence[Stream]:
        """List of all inlet streams."""
        return self._ins    
    @property
    def outs(self) -> Sequence[Stream]:
        """List of all outlet streams."""
        return self._outs
    
    @property
    def feed(self) -> Stream:
        """Equivalent to ins[0] when the number of inlets is 1."""
        streams = self._ins._streams
        size = len(streams)
        if size == 1: return streams[0]
        elif size > 1: raise AttributeError(f"{repr(self)} has more than one inlet")
        else: raise AttributeError(f"{repr(self)} has no inlet")
    @feed.setter
    def feed(self, feed): 
        ins = self._ins
        streams = ins._streams
        size = len(streams)
        if size == 1: ins[0] = feed
        elif size > 1: raise AttributeError(f"{repr(self)} has more than one inlet")
        else: raise AttributeError(f"{repr(self)} has no inlet")
    inlet = influent = feed
    
    @property
    def product(self) -> Stream:
        """Equivalent to outs[0] when the number of outlets is 1."""
        streams = self._outs._streams
        size = len(streams)
        if size == 1: return streams[0]
        elif size > 1: raise AttributeError(f"{repr(self)} has more than one outlet")
        else: raise AttributeError(f"{repr(self)} has no outlet")
    @product.setter
    def product(self, product): 
        outs = self._outs
        streams = outs._streams
        size = len(streams)
        if size == 1: outs[0] = product
        elif size > 1: raise AttributeError(f"{repr(self)} has more than one outlet")
        else: raise AttributeError(f"{repr(self)} has no outlet")
    outlet = effluent = product
    
    @property
    def flat_ins(self) -> Sequence[Stream]:
        """List of all inlet single-phase streams."""
        flat = []
        for i in self._ins:
            if len(i.phases) > 1: flat.extend(i)
            else: flat.append(i)
        return flat
    @property
    def flat_outs(self) -> Sequence[Stream]:
        """List of all outlet single-phase streams."""
        flat = []
        for i in self._outs:
            if len(i.phases) > 1: flat.extend(i)
            else: flat.append(i)
        return flat
        
    @property
    def auxiliary_units(self) -> list[Unit]:
        """Return list of all auxiliary units."""
        getfield = getattr
        isa = isinstance
        auxiliary_units = []
        for name in self.auxiliary_unit_names:
            unit = getfield(self, name, None)
            if unit is None: continue 
            if isa(unit, Iterable):
                auxiliary_units.extend(unit)
            else:
                auxiliary_units.append(unit)
        return auxiliary_units

    @property
    def nested_auxiliary_units(self) -> list[Unit]:
        """Return list of all auxiliary units, including nested ones."""
        getfield = getattr
        isa = isinstance
        auxiliary_units = []
        for name in self.auxiliary_unit_names:
            unit = getfield(self, name, None)
            if unit is None: continue 
            if isa(unit, Iterable):
                auxiliary_units.extend(unit)
                for u in unit:
                    if not isinstance(u, Unit): continue
                    for auxunit in u.auxiliary_units:
                        auxiliary_units.append(auxunit)
                        auxiliary_units.extend(auxunit.nested_auxiliary_units)
            else:
                auxiliary_units.append(unit)
                if not isinstance(unit, Unit): continue
                for auxunit in unit.auxiliary_units:
                    auxiliary_units.append(auxunit)
                    auxiliary_units.extend(auxunit.nested_auxiliary_units)
        return auxiliary_units

    def _diagram_auxiliary_units_with_names(self) -> list[tuple[str, Unit]]:
        """Return list of name - auxiliary unit pairs."""
        getfield = getattr
        isa = isinstance
        auxiliary_units = []
        names = (
            self.diagram_auxiliary_unit_names 
            if hasattr(self, 'diagram_auxiliary_unit_names')
            else self.auxiliary_unit_names
        )
        for name in names:
            unit = getfield(self, name, None)
            if unit is None: continue 
            if isa(unit, Iterable):
                for i, u in enumerate(unit):
                    auxiliary_units.append(
                        (f"{name}[{i}]", u)
                    )
            else:
                auxiliary_units.append(
                    (name, unit)
                )
        return auxiliary_units

    def get_auxiliary_units_with_names(self) -> list[tuple[str, Unit]]:
        """Return list of name - auxiliary unit pairs."""
        getfield = getattr
        isa = isinstance
        auxiliary_units = []
        for name in self.auxiliary_unit_names:
            unit = getfield(self, name, None)
            if unit is None: continue 
            if isa(unit, Iterable):
                for i, u in enumerate(unit):
                    auxiliary_units.append(
                        (f"{name}[{i}]", u)
                    )
            else:
                auxiliary_units.append(
                    (name, unit)
                )
        return auxiliary_units

    def _diagram_nested_auxiliary_units_with_names(self, depth=-1) -> list[Unit]:
        """Return list of all diagram auxiliary units, including nested ones."""
        auxiliary_units = []
        if depth: 
            depth -= 1
        else:
            return auxiliary_units
        for name, auxunit in self._diagram_auxiliary_units_with_names():
            if auxunit is None: continue 
            auxiliary_units.append((name, auxunit))
            if not isinstance(auxunit, Unit): continue
            auxiliary_units.extend(
                [('.'.join([name, i]), j)
                 for i, j in auxunit._diagram_nested_auxiliary_units_with_names(depth)]
            )
        return auxiliary_units

    def get_nested_auxiliary_units_with_names(self, depth=-1) -> list[Unit]:
        """Return list of all auxiliary units, including nested ones."""
        auxiliary_units = []
        if depth: 
            depth -= 1
        else:
            return auxiliary_units
        for name, auxunit in self.get_auxiliary_units_with_names():
            if auxunit is None: continue 
            auxiliary_units.append((name, auxunit))
            if not isinstance(auxunit, Unit): continue
            auxiliary_units.extend(
                [('.'.join([name, i]), j)
                 for i, j in auxunit.get_nested_auxiliary_units_with_names(depth)]
            )
        return auxiliary_units

    def _unit_auxins(self, N_streams, streams, thermo):
        if streams is None or streams == ():
            return [self.auxin(Stream(None, thermo=thermo)) for i in range(N_streams)]
        elif isinstance(streams, (Stream, MissingStream)) or streams.__class__ is str:
            return self.auxin(streams, thermo=thermo)
        else:
            return [self.auxin(i, thermo=thermo) for i in streams]

    def _unit_auxouts(self, N_streams, streams, thermo):
        if streams is None or streams == ():
            return [self.auxout(Stream(None, thermo=thermo)) for i in range(N_streams)]
        elif isinstance(streams, (Stream, MissingStream)) or streams.__class__ is str:
            return self.auxout(streams, thermo=thermo)
        else:
            return [self.auxout(i, thermo=thermo) for i in streams]
    
    def register_auxiliary(self, unit, name=None):
        names_list = self.auxiliary_unit_names
        names_set = set(names_list)
        if isinstance(unit, Unit):
            if name is None: name = unit._ID
            if name not in names_set: names_list.append(name)
            setattr(self, name, unit)
        elif isinstance(unit, Iterable):
            if name is None: raise ValueError('`name` must be a string')
            if name not in names_set: names_list.append(name)
            setattr(self, name, unit)
        elif hasattr(unit, 'units'):
            system = unit
            if name is None: name = system._ID
            setattr(self, name, system)
            for i in system.units:
                name = i._ID
                if name not in names_set: names_list.append(name)
                setattr(self, name, i)
        else:
            raise ValueError('`unit` must be a unit, list of units, or a system')
    
    def auxiliary(
            self, name, cls, ins=None, outs=(), thermo=None,
            register=False, **kwargs
        ):
        """
        Create and register an auxiliary unit operation. Inlet and outlet
        streams automatically become auxlets so that parent unit streams will
        not disconnect.

        """
        if thermo is None: thermo = self.thermo
        auxunit = cls.__new__(cls)
        stack = getattr(self, name, None)
        if isinstance(stack, list): 
            name = f"{name}[{len(stack)}]"
            stack.append(auxunit)
        else:
            setattr(self, name, auxunit)
        if register and name not in self.auxiliary_unit_names:
           self.auxiliary_unit_names.add(name)
        auxunit.owner = self # Avoids property package checks
        auxunit.__init__(
            '.' + name, 
            self._unit_auxins(cls._N_ins, ins, thermo), 
            self._unit_auxouts(cls._N_outs, outs, thermo),
            thermo, 
            **kwargs
        )
        return auxunit

    def auxlet(self, stream: Stream, thermo=None):
        """
        Define auxiliary unit inlet or outlet. This method has two
        behaviors:

        * If the stream is not connected to this unit, define the Stream 
          object's source or sink to be this unit without actually connecting 
          it to this unit.

        * If the stream is already connected to this unit, return a superposition
          stream which can be connected to auxiliary units without being disconnected
          from this unit. 

        """
        if thermo is None: thermo = self.thermo
        if stream is None: 
            stream = Stream(None, thermo=thermo)
            stream._sink = stream._source = self
        elif isinstance(stream, str): 
            stream = Stream('.' + stream, thermo=thermo)
            stream._source = stream._sink = self
        elif self is stream._source and stream in self._outs:
            port = OutletPort.from_outlet(stream)
            stream = SuperpositionOutlet(port)
            self.auxouts[port.index] = stream
        elif self is stream._sink and stream in self._ins:
            port = InletPort.from_inlet(stream)
            stream = SuperpositionInlet(port)
            self.auxins[port.index] = stream
        else:
            if stream._source is None: stream._source = self
            if stream._sink is None: stream._sink = self
        return stream

    def auxin(self, stream: Stream, thermo=None):
        """
        Define auxiliary unit inlet. This method has two
        behaviors:

        * If the stream is not connected to this unit, define the Stream 
          object's source to be this unit without actually connecting it to this unit.

        * If the stream is already connected to this unit, return a superposition
          stream which can be connected to auxiliary units without being disconnected
          from this unit.

        """
        if thermo is None: thermo = self.thermo
        if stream is None: 
            stream = Stream(None, thermo=thermo)
            stream._sink = stream._source = self
        elif isinstance(stream, str): 
            stream = Stream('.' + stream, thermo=thermo)
            stream._source = stream._sink = self
        elif self is stream._sink and stream in self._ins:
            port = InletPort.from_inlet(stream)
            stream = SuperpositionInlet(port)
            self.auxins[port.index] = stream
        else:
            if stream._source is None: stream._source = self
            if stream._sink is None: stream._sink = self
        return stream

    def auxout(self, stream: Stream, thermo=None):
        """
        Define auxiliary unit outlet. This method has two
        behaviors:

        * If the stream is not connected to this unit, define the Stream 
          object's source or sink to be this unit without actually connecting 
          it to this unit.

        * If the stream is already connected to this unit, return a superposition
          stream which can be connected to auxiliary units without being disconnected
          from this unit. 

        """
        if thermo is None: thermo = self.thermo
        if stream is None: 
            stream = Stream(None, thermo=thermo)
            stream._sink = stream._source = self
        elif isinstance(stream, str): 
            stream = Stream('.' + stream, thermo=thermo)
            stream._source = stream._sink = self
        elif self is stream._source and stream in self._outs:
            port = OutletPort.from_outlet(stream)
            stream = SuperpositionOutlet(port)
            self.auxouts[port.index] = stream
        else:
            if stream._source is None: stream._source = self
            if stream._sink is None: stream._sink = self
        return stream

    def _assembled_from_auxiliary_units(self):
        #: Serves for checking whether to include this unit in auxiliary diagrams.
        #: If all streams are in common, it must be assembled by auxiliary units.
        return not set([i.ID for i in self.ins + self.outs]).difference(
            sum([[i.ID for i in (i.ins + i.outs)] for i in self.auxiliary_units], [])
        )
    
    def get_node(self):
        """Return unit node attributes for graphviz."""
        if tmo.preferences.minimal_nodes:
            return self._graphics.get_minimal_node(self)
        else:
            node = self._graphics.get_node_tailored_to_unit(self)
            try:
                node['tooltip'] = self._get_tooltip_string()
            except:
                pass
            return node
    
    def add_specification(self, 
            f: Optional[Callable]=None, 
            run: Optional[bool]=None, 
            args: Optional[tuple]=(),
            impacted_units: Optional[tuple[Unit, ...]]=None,
            prioritize: Optional[bool]=None,
        ):
        """
        Add a specification.

        Parameters
        ----------
        f : 
            Specification function runned for mass and energy balance.
        run : 
            Whether to run the built-in mass and energy balance after 
            specifications. Defaults to False.
        args : 
            Arguments to pass to the specification function.
        impacted_units :
            Other units impacted by specification. The system will make sure to 
            run itermediate upstream units when simulating.
        prioritize :
            Whether to prioritize the unit operation within a recycle loop (if any).
            
        Examples
        --------
        :doc:`../tutorial/Process_specifications`

        See Also
        --------
        add_bounded_numerical_specification
        specifications
        run

        Notes
        -----
        This method also works as a decorator.

        """
        if not f: return lambda f: self.add_specification(f, run, args, impacted_units, prioritize)
        if not callable(f): raise ValueError('specification must be callable')
        self._specifications.append(ProcessSpecification(f, args, impacted_units))
        if run is not None: self.run_after_specifications = run
        if prioritize is not None: self.prioritize = prioritize
        return f
    
    def add_bounded_numerical_specification(self, f=None, *args, **kwargs):
        """
        Add a bounded numerical specification that solves x where f(x) = 0 using an 
        inverse quadratic interpolation solver.
        
        Parameters
        ----------
        f : Callable
            Objective function in the form of f(x, *args).
        x : float, optional
            Root guess.
        x0, x1 : float
            Root bracket. Solution must lie within x0 and x1.
        xtol : float, optional
            Solver stops when the root lies within xtol. Defaults to 0.
        ytol : float, optional 
            Solver stops when the f(x) lies within ytol of the root. Defaults to 5e-8.
        args : tuple, optional
            Arguments to pass to f.
        maxiter : 
            Maximum number of iterations. Defaults to 50.
        checkiter : bool, optional
            Whether to raise a Runtime error when tolerance could not be 
            satisfied before the maximum number of iterations. Defaults to True.
        checkroot : bool, optional
            Whether satisfying both tolerances, xtol and ytol, are required 
            for termination. Defaults to False.
        checkbounds : bool, optional
            Whether to raise a ValueError when in a bounded solver when the 
            root is not certain to lie within bounds (i.e. f(x0) * f(x1) > 0.).
            Defaults to True.
            
        Examples
        --------
        :doc:`../tutorial/Process_specifications`

        See Also
        --------
        add_specification
        specifications
        
        Notes
        -----
        This method also works as a decorator.

        """
        if not f: return lambda f: self.add_bounded_numerical_specification(f, *args, **kwargs)
        if not callable(f): raise ValueError('f must be callable')
        self._specifications.append(BoundedNumericalSpecification(f, *args, **kwargs))
        return f
    
    def run(self):
        """
        Run mass and energy balance. This method also runs specifications
        user defined specifications unless it is being run within a 
        specification (to avoid infinite loops). 
        
        See Also
        --------
        _run
        specifications
        add_specification
        add_bounded_numerical_specification
        
        """
        if self._skip_simulation_when_inlets_are_empty and all([i.isempty() for i in self._ins]):
            for i in self._outs: i.empty()
            return
        self._run_with_specifications()
        
    def _run_with_specifications(self):
        specifications = self._specifications
        if specifications:
            active_specifications = self._active_specifications
            if len(active_specifications) == len(specifications):
                self._run()
            else:
                for ps in specifications: 
                    if ps in active_specifications: continue
                    active_specifications.add(ps)
                    try: ps()
                    finally: active_specifications.remove(ps)
                if self.run_after_specifications: self._run()
        else:
            self._run()
    
    def path_from(self, units, inclusive=False, system=None):
        """
        Return a tuple of units and systems starting from `units` until this one 
        (not inclusive by default).
        
        """
        units = (units,) if isinstance(units, Unit) else tuple(units)
        if system: # Includes recycle loops
            path = system.path_section(units, (self,))
        else: # Path outside system, so recycle loops may not converge (and don't have to)
            path = []
            added_units = set()
            upstream_units = self.get_upstream_units()
            for unit in units:
                if unit in upstream_units: add_path_segment(unit, self, path, added_units)
            if inclusive and unit not in added_units: path.append(unit)
        return path        
    
    def path_until(self, units, inclusive=False, system=None):
        """
        Return a tuple of units and systems starting from this one until the end
        units (not inclusive by default).
        
        """
        units = (units,) if isinstance(units, Unit) else tuple(units)
        if system: # Includes recycle loops
            path = system.path_section((self,), units, inclusive)
        else: # Path outside system, so recycle loops may not converge (and don't have to)
            path = []
            added_units = set()
            downstream_units = self.get_downstream_units()
            for unit in units:
                if unit in downstream_units: add_path_segment(self, unit, path, added_units)
            if inclusive and unit not in added_units: path.append(unit)
        return path
    
    def run_until(self, units, inclusive=False, system=None):
        """
        Run all units and converge all systems starting from this one until the end units
        (not inclusive by default).
        
        See Also
        --------
        path_until
        
        """
        isa = isinstance
        for i in self.path_until(units, inclusive, system): 
            if isa(i, Unit): i.run()
            else: i.converge() # Must be a system
        
    @property
    def system(self):
        return self._system
    
    @ignore_docking_warnings
    def disconnect(self, discard=False, inlets=None, outlets=None, join_ends=False):
        ins = self._ins
        outs = self._outs
        if inlets is None: 
            inlets = ins._streams.copy()
            ins[:] = ()
        else:
            for i in inlets: ins[ins.index(i) if isinstance(i, Stream) else i] = None
        if outlets is None: 
            outlets = outs._streams.copy()
            outs[:] = ()
        else:
           for o in outlets: outs[outs.index(o) if isinstance(o, Stream) else o] = None
        if join_ends:
            if len(inlets) != len(outlets):
                inlets = [i for i in inlets if i]
                outlets = [i for i in outlets if i]
                if len(inlets) != len(outlets):
                    raise ValueError("number of inlets must match number of outlets to join ends")
            for inlet, outlet in zip(inlets, outlets):
                if outlet.sink: outlet.sink.ins.replace(outlet, inlet)
        if discard: self.registry.discard(self)

    @ignore_docking_warnings
    def insert(self, stream: Stream, inlet: int|Stream=None, outlet: int|Stream=None):
        """
        Insert unit between two units at a given stream connection.
        
        Examples
        --------
        >>> from biosteam import *
        >>> settings.set_thermo(['Water'], cache=True)
        >>> feed = Stream('feed')
        >>> other_feed = Stream('other_feed')
        >>> P1 = Pump('P1', feed, 'pump_outlet')
        >>> H1 = HXutility('H1', P1-0, T=310)
        >>> M1 = Mixer('M1', other_feed, 'mixer_outlet')
        >>> M1.insert(P1-0)
        >>> M1.show()
        Mixer: M1
        ins...
        [0] other_feed
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        [1] pump_outlet  from  Pump-P1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        outs...
        [0] mixer_outlet  to  HXutility-H1
            phase: 'l', T: 298.15 K, P: 101325 Pa
            flow: 0
        
        """
        source = stream.source
        sink = stream.sink
        added_unit = False
        if outlet is None:
            if self._outs_size_is_fixed:
                if self._N_outs == 1:
                    sink.ins.replace(stream, self.outs[0])
                else:
                    raise ValueError("undefined outlet; must pass outlet when outlets are fixed and multiple are available")
            else:
                self.outs.append(stream)
                added_unit = True
        else:
            if isinstance(outlet, Stream):
                if outlet.source is not self:
                    raise ValueError("source of given outlet must be this object")
            else:
                outlet = self.outs[outlet]
            sink.ins.replace(stream, outlet)
        if inlet is None:
            if self._ins_size_is_fixed or added_unit:
                if self._N_ins == 1:
                    source.outs.replace(stream, self.ins[0])
                else:
                    raise ValueError("undefined inlet; must pass inlet when inlets are fixed and multiple are available")
            else:
                self.ins.append(stream)
        else:
            if isinstance(inlet, Stream):
                if inlet.sink is not self:
                    raise ValueError("sink of given inlet must be this object")
            else:
                inlet = self.outs[inlet]
            source.outs.replace(stream, inlet)
    
    @ignore_docking_warnings
    def take_place_of(self, other, discard=False):
        """Replace inlets and outlets from this unit with that of 
        another unit."""
        self._ins[:] = other.ins
        self._outs[:] = other.outs
        if discard: self.registry.discard(other)
    
    @ignore_docking_warnings
    def replace_with(self, other=None, discard=False):
        """Replace inlets and outlets from another unit with this unit."""
        if other is None:
            ins = self._ins
            outs = self._outs
            for inlet, outlet in zip(tuple(ins), tuple(outs)):
                source = inlet.source
                if source:
                    source.outs.replace(inlet, outlet)
                else:
                    sink = outlet.sink
                    if sink: sink.ins.replace(outlet, inlet)
            ins.empty()
            outs.empty()
        else:
            other.ins[:] = self._ins
            other.outs[:] = self._outs
        if discard: self.registry.discard(self)
                    
    # Forward pipping
    def __sub__(self, other):
        """Source streams."""
        isa = isinstance
        if hasattr(other, '_ins'):
            other._ins[:] = self._outs
            return other
        elif isa(other, int_types):
            return self._outs[other]
        elif isa(other, Stream):
            self._outs[:] = (other,)
            return self
        elif isa(other, (tuple, list, np.ndarray)):
            if all([isa(i, int_types) for i in other]):
                outs = self._outs
                return [outs[i] for i in other]
            else:
                self._outs[:] = other
                return self
        else:
            return other.__rsub__(self)

    def __rsub__(self, other):
        """Sink streams."""
        isa = isinstance
        if isa(other, int_types):
            return self._ins[other]
        elif isa(other, Stream):
            self._ins[:] = (other,)
            return self
        elif isa(other, (tuple, list, np.ndarray)):
            if all([isa(i, int_types) for i in other]):
                ins = self._ins
                return [ins[i] for i in other]
            else:
                self._ins[:] = other
                return self
        else:
            raise ValueError(f"cannot pipe '{type(other).__name__}' object")

    
    def _add_upstream_neighbors_to_set(self, set, ends, universal):
        """Add upsteam neighboring units to set."""
        for s in self._ins:
            u = s._source
            if u and (universal or not u._universal) and not (ends and s in ends):
                set.add(u)

    def _add_downstream_neighbors_to_set(self, set, ends, universal):
        """Add downstream neighboring units to set."""
        for s in self._outs:
            u = s._sink
            if u and (universal or not u._universal) and not (ends and s in ends):
                set.add(u)

    def get_downstream_units(self, ends=None, universal=True, downstream_units=None):
        """Return a set of all units downstream."""
        if downstream_units is None: 
            downstream_units = set()
            new_length = 0
        else:
            new_length = len(downstream_units)
        outer_periphery = set()
        self._add_downstream_neighbors_to_set(outer_periphery, ends, universal)
        inner_periphery = None
        old_length = -1
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            downstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                unit._add_downstream_neighbors_to_set(outer_periphery, ends, universal)
            new_length = len(downstream_units)
        return downstream_units
    
    def get_upstream_units(self, ends=None, universal=True, upstream_units=None):
        """Return a set of all units upstream."""
        if upstream_units is None:
            upstream_units = set()
            new_length = 0
        else:
            new_length = len(upstream_units)
        outer_periphery = set()
        self._add_upstream_neighbors_to_set(outer_periphery, ends, universal)
        inner_periphery = None
        old_length = -1
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            upstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                unit._add_upstream_neighbors_to_set(outer_periphery, ends, universal)
            new_length = len(upstream_units)
        return upstream_units
    
    def get_recycle_units(self):
        downstream_units = set()
        downstream_ends = set()
        upstream_units = set()
        upstream_ends = set()
        self.get_downstream_units(
            downstream_units=downstream_units, 
            ends=downstream_ends, 
            universal=False,
        )
        self.get_upstream_units(
            upstream_units=upstream_units, 
            ends=upstream_ends, 
            universal=False,
        )
        return downstream_units.intersection(upstream_units)
    
    def neighborhood(self, 
            radius: Optional[int]=1, 
            upstream: Optional[bool]=True,
            downstream: Optional[bool]=True, 
            ends: Optional[Stream]=None, 
            universal: Optional[bool]=None
        ):
        """
        Return a set of all neighboring units within given radius.
        
        Parameters
        ----------
        radius : 
            Maximum number streams between neighbors.
        downstream : 
            Whether to include downstream units.
        upstream : 
            Whether to include upstream units.
        ends :
            Streams that mark the end of the neighborhood.
        universal :
            Whether to include universal units.
        
        """
        radius -= 1
        neighborhood = set()
        if radius < 0: return neighborhood
        if upstream:self._add_upstream_neighbors_to_set(neighborhood, ends, universal)
        if downstream: self._add_downstream_neighbors_to_set(neighborhood, ends, universal)
        direct_neighborhood = neighborhood
        for i in range(radius):
            neighbors = set()
            for unit in direct_neighborhood:
                if upstream: unit._add_upstream_neighbors_to_set(neighbors, ends, universal)
                if downstream: unit._add_downstream_neighbors_to_set(neighbors, ends, universal)
            if neighbors == direct_neighborhood: break
            direct_neighborhood = neighbors
            neighborhood.update(direct_neighborhood)
        return neighborhood

    def diagram(self, radius: Optional[int]=0, upstream: Optional[bool]=True,
                downstream: Optional[bool]=True, file: Optional[str]=None, 
                format: Optional[str]=None, display: Optional[bool]=True,
                auxiliaries: Optional[bool]=-1,
                **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram
        of the unit and all neighboring units within given radius.
        
        Parameters
        ----------
        radius : 
            Maximum number streams between neighbors.
        downstream : 
            Whether to show downstream operations.
        upstream : 
            Whether to show upstream operations.
        file : 
            Must be one of the following:
            
            * [str] File name to save diagram.
            * [None] Display diagram in console.
            
        format : 
            Format of file.
        display : 
            Whether to display diagram in console or to return the graphviz 
            object.
        auxiliaries:
            Depth of auxiliary units to display.
        
        """
        if radius > 0:
            units = self.neighborhood(radius, upstream, downstream)
            units.add(self)
        else:
            units = [self]
        return bst.System(None, units).diagram(
            format=format, auxiliaries=auxiliaries, display=display, 
            file=file, title='', **graph_attrs
        )
    
    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    # Representation
    def _info(self, layout, T, P, flow, composition, N, IDs, sort, data):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        ins = [(i.materialize() if hasattr(i, 'materialize') else i) for i in self.ins]
        outs = [(i.materialize() if hasattr(i, 'materialize') else i) for i in self.outs]
        return info + repr_ins_and_outs(layout, ins, outs, T, P, flow, composition, N, IDs, sort, data)

    def show(self, layout=None, T=None, P=None, flow=None, composition=None, N=None, IDs=None, sort=None, data=True):
        """Prints information on unit."""
        print(self._info(layout, T, P, flow, composition, N, IDs, sort, data))

    def _ipython_display_(self):
        if bst.preferences.autodisplay: self.diagram()
        self.show()

    def _init_ins(self, ins):
        self.auxins = {}
        self._init_inlets(ins)
    
    def _init_outs(self, outs):
        self.auxouts = {}
        self._init_outlets(outs)

    def _init_utils(self):
        self.heat_utilities = [HeatUtility for i in range(getattr(self, '_N_heat_utilities', 0))]
        self.power_utility = PowerUtility()
        
    def _init_results(self):
        try: self.F_BM = self._F_BM_default.copy()
        except AttributeError: self.F_BM = {}
        self.F_D = {}
        self.F_P = {}
        self.F_M = {}
        self.design_results = {}
        self.baseline_purchase_costs = {}
        self.purchase_costs = {}
        self.installed_costs = {}
        self._inlet_utility_indices = {}
        self._outlet_utility_indices = {}
        self._inlet_cost_indices = {}
        self._outlet_revenue_indices = {}
        self._flow_fractions = {}
        try: self.equipment_lifetime = copy(self._default_equipment_lifetime)
        except AttributeError: self.equipment_lifetime = {}

    def _init_specifications(self):
        self._specifications = []
        self.run_after_specifications = False
        self._active_specifications = set()
    
    def _reset_thermo(self, thermo):
        for i in (self._ins._streams + self._outs._streams):
            try:
                if i: i._reset_thermo(thermo)
            except:
                raise RuntimeError(f'failed to reset {repr(self)}.thermo')
        if thermo is self.thermo: return
        self._load_thermo(thermo)
        chemicals = thermo.chemicals
        isa = isinstance
        hasfield = hasattr
        def reset_thermo(obj, old=set()):
            hash = id(obj)
            if hash in old: return 
            old.add(hash)
            if isa(obj, tmo.ReactionSystem):
                for rxn in obj._reactions:
                    if hasfield(rxn, 'reset_chemicals') and rxn.chemicals is not chemicals:
                        rxn.reset_chemicals(chemicals)
            elif hasfield(obj, 'reset_chemicals') and obj.chemicals is not chemicals:
                obj.reset_chemicals(chemicals)
            elif hasfield(obj, '_reset_thermo') and obj.thermo is not thermo:
                obj._reset_thermo(thermo)
            elif isa(obj, dict):
                for i in obj.values(): reset_thermo(i)
            elif isa(obj, Iterable):
                for i in obj: reset_thermo(i)
                
        for obj in self.__dict__.values(): reset_thermo(obj)
            
    def response(self, name):
        """Register response for convergence model prediction."""
        self.responses.add(
            bst.GenericResponse(self, name)
        )
            
    @property
    def net_power(self) -> float:
        """Net power consumption [kW]."""
        return self.power_utility.power
    @property
    def net_duty(self) -> float:
        """Net duty including heat transfer losses [kJ/hr]."""
        return sum([i.duty for i in self.heat_utilities])
    @property
    def net_cooling_duty(self) -> float:
        """Net cooling duty including heat transfer losses [kJ/hr]."""
        return sum([i.duty for i in self.heat_utilities if i.duty < 0.])
    @property
    def net_heating_duty(self) -> float:
        """Net cooling duty including heat transfer losses [kJ/hr]."""
        return sum([i.duty for i in self.heat_utilities if i.duty > 0.])
    
    def add_power_utility(self, power):
        """Add power utility [kW]. Use a positive value for consumption and 
        a negative for production."""
        power_utility = self.power_utility
        if power >= 0.:
            power_utility.consumption += power
        else:
            power_utility.production -= power
    
    def create_heat_utility(self,
            agent: Optional[UtilityAgent]=None,
            heat_transfer_efficiency: Optional[float]=None,
        ):
        """Create heat utility object associated to unit."""
        hu = HeatUtility(heat_transfer_efficiency, None)
        self.heat_utilities.append(hu)
        return hu
    
    def add_heat_utility(self, 
            unit_duty: float, 
            T_in: float,
            T_out: Optional[float]=None, 
            agent: Optional[UtilityAgent]=None,
            heat_transfer_efficiency: Optional[float]=None,
            hxn_ok: Optional[bool]=False,
        ):
        """
        Add utility requirement given the duty and inlet and outlet 
        temperatures.
        
        Parameters
        ----------
        unit_duty :
               Unit duty requirement [kJ/hr]
        T_in : 
               Inlet process stream temperature [K]
        T_out : 
               Outlet process stream temperature [K]
        agent : 
                Utility agent to use. Defaults to a suitable agent from 
                predefined heating/cooling utility agents.
        heat_transfer_efficiency : 
            Enforced fraction of heat transferred from utility (due
            to losses to environment).
        hxn_ok :
            Whether heat utility can be satisfied within a heat exchanger network.
            
        """
        hu = HeatUtility(heat_transfer_efficiency, self, hxn_ok)
        self.heat_utilities.append(hu)
        hu(unit_duty, T_in, T_out, agent)
        return hu
    
    def define_utility(self, name: str, stream: Stream):
        """
        Define an inlet or outlet stream as a utility by name.
        
        Parameters
        ----------
        name : 
            Name of utility, as defined in :meth:`settings.stream_prices <thermosteam._settings.ProcessSettings.stream_prices>`.
        stream :
            Inlet or outlet utility stream.
        
        """
        if name not in bst.stream_prices:
            raise ValueError(f"price of '{name}' must be defined in settings.stream_prices")
        if stream._sink is self:
            self._inlet_utility_indices[name] = self._ins._streams.index(stream)
        elif stream._source is self:
            self._outlet_utility_indices[name] = self._outs._streams.index(stream)
        else:
            raise ValueError(f"stream '{stream.ID}' must be connected to {repr(self)}")
    
    def define_credit(self, name: str, stream: Stream, fraction: float|Callable=None):
        """
        Define an inlet or outlet stream as a fee/credit by name.
        
        Parameters
        ----------
        name : 
            Name of fee/credit, as defined in :meth:`settings.stream_prices <thermosteam._settings.ProcessSettings.stream_prices>`.
        stream :
            Inlet or outlet fee/credit stream.
        fraction :
            Fraction of the stream that receives the fee/credit. If given a function,
            it must return the fraction.
        
        """
        if name not in bst.stream_prices:
            raise ValueError(f"price of '{name}' must be defined in settings.stream_prices")
        if stream._sink is self:
            self._inlet_cost_indices[name] = self._ins._streams.index(stream)
        elif stream._source is self:
            self._outlet_revenue_indices[name] = self._outs._streams.index(stream)
        else:
            raise ValueError(f"stream '{stream.ID}' must be connected to {repr(self)}")
        if fraction: self._flow_fractions[name] = fraction
    
    define_fee = define_credit
    
    def _get_flow_fractions(self):
        """Return fee/credit fractions with any callables evaluated."""
        return {name: fraction() if callable(fraction) else fraction
                for name, fraction in self._flow_fractions.items()}

    def get_inlet_cost_flows(self):
        ins = self._ins._streams
        flows = {name: ins[index].F_mass for name, index in (self._inlet_utility_indices | self._inlet_cost_indices).items()}
        for name, fraction in self._flow_fractions.items(): flows[name] *= fraction() if callable(fraction) else fraction
        return flows

    def get_outlet_revenue_flows(self):
        outs = self._outs._streams
        flows = {name: outs[index].F_mass for name, index in (self._outlet_utility_indices | self._outlet_revenue_indices).items()}
        for name, fraction in self._flow_fractions.items(): flows[name] *= fraction() if callable(fraction) else fraction
        return flows
    
    def get_design_and_capital(self):
        return UnitDesignAndCapital(
            self, self.F_BM.copy(), self.F_D.copy(), self.F_P.copy(), self.F_M.copy(), 
            self.design_results.copy(), self.baseline_purchase_costs.copy(),
            self.purchase_costs.copy(), self.installed_costs.copy(),
        )
    
    def get_agile_design_and_capital(self, design_and_capital: list[UnitDesignAndCapital]):
        names = (
            'F_BM', 'F_D', 'F_P', 'F_M', 'design_results',
            'baseline_purchase_costs', 'purchase_costs', 'installed_costs',
        )
        max_agile_design = getattr(self, '_max_agile_design', None)
        agile_scenario = {i: {} for i in names}
        if max_agile_design:
            designs = [i.design_results for i in design_and_capital]
            agile_design = agile_scenario['design_results']
            for design_results in designs:
                for i in max_agile_design:
                    if i in design_results:
                        j = design_results[i]
                    else:
                        continue
                    if i in agile_design:
                        if abs(j) > abs(agile_design[i]):
                            agile_design[i] = j
                    else:
                        agile_design[i] = j
            names = ('design_results', 'baseline_purchase_costs', 
                     'purchase_costs', 'installed_costs')
            dcts = [getattr(self, i) for i in names]
            self.design_results = agile_design
            for i in names[1:]: setattr(self, i, {}) 
            Unit._setup(self)
            try:
                self._cost()
                self._load_costs()
            except:
                warn(f"failed to create agile design for {repr(self)}; "
                      "assuming design with highest capital cost will do",
                      category=RuntimeWarning, stacklevel=2)
            else:
                design_and_capital.append(self.get_design_and_capital())
            finally:
                for i, j in zip(names, dcts): setattr(self, i, j) 
        self._fill_agile_design_and_capital(agile_scenario, design_and_capital)
        agile_scenario['unit'] = self
        return UnitDesignAndCapital(**agile_scenario)
                
    def _fill_agile_design_and_capital(self, agile_scenario, design_and_capital):
        for results in design_and_capital:
            for name, maxdct in agile_scenario.items():
                if name == 'design_results': continue
                dct = getattr(results, name)
                for i, j in dct.items():
                    if i in maxdct:
                        if abs(j) > abs(maxdct[i]):
                            maxdct[i] = j
                    else:
                        maxdct[i] = j
        F_BM = agile_scenario['F_BM']
        F_D = agile_scenario['F_D']
        F_P = agile_scenario['F_P']
        F_M = agile_scenario['F_M']
        baseline_purchase_costs = agile_scenario['baseline_purchase_costs']
        purchase_costs = agile_scenario['purchase_costs']
        installed_costs = agile_scenario['installed_costs']
        materials_and_maintenance = self._materials_and_maintenance
        for name, Cp in baseline_purchase_costs.items():
            F = F_D.get(name, 1.) * F_P.get(name, 1.) * F_M.get(name, 1.)
            installed_cost = Cp * (F_BM.get(name, 1.) + F - 1.)
            purchase_cost = Cp * F
            if installed_cost > installed_costs[name]:
                if name in materials_and_maintenance:
                    installed_costs[name] += installed_cost
                else:
                    installed_costs[name] = installed_cost
            if purchase_cost > purchase_costs[name]:
                if name in materials_and_maintenance:
                    purchase_costs[name] += purchase_cost
                else:
                    purchase_costs[name] = purchase_cost
    
    def _assert_compatible_property_package(self):
        if self.owner is not self: return
        CASs = self.chemicals.CASs
        streams = self._ins + self._outs
        assert all([s.chemicals.CASs == CASs for s in streams if s]), (
            "unit operation chemicals are incompatible with inlet and outlet streams; "
            "try using the `thermo` keyword argument to initialize the unit operation "
            "with a compatible thermodynamic property package"
        )
    
    def _load_costs(self):
        r"""
        Calculate and save free on board (f.o.b.) purchase costs and
        installed equipment costs (i.e. bare-module cost) for each item in the 
        :attr:`~Unit.baseline_purchase_costs` dictionary and in auxiliary units. This 
        method is run after the :attr:`~Unit._cost` method at the end of unit simulation.
        
        Notes
        -----
        As explained in [1]_, the f.o.b. purchase cost is given by:
        
        .. math::
           
           C_{P} = C_{Pb}F_{D}F_{P}F_{M}
        
        And the installed equipment cost is given by:
        
        .. math::
           
           C_{BM} = C_{Pb} (F_{BM} + F_{D}F_{P}F_{M} - 1)
        
        Where:
            * :math:`C_{Pb}`: Baseline purchase cost.
            * :math:`F_{BM}`: Bare module factor.
            * :math:`F_{D}`: Design factor.
            * :math:`F_{P}`: Pressure factor.
            * :math:`F_{M}`: Material factor.
        
        Values for the bare-module, design, pressure, and material factors of 
        each equipment should be stored in the :attr:`~Unit.F_BM`, :attr:`~Unit.F_D`, 
        :attr:`~Unit.F_P`, and :attr:`~Unit.F_M` dictionaries.
        
        Warning
        -------
        If an item is listed in the :attr:`~Unit.purchase_costs` dictionary but not in the
        :attr:`~Unit.baseline_purchase_costs` dictionary, the baseline purchase cost is 
        assumed to be the same as the purchase cost.
        
        References
        ----------
        .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley. Cost Accounting and Capital Cost Estimation (Chapter 16)
        
        """
        if getattr(self, '_costs_loaded', False): return
        F_BM = self.F_BM
        F_D = self.F_D
        F_P = self.F_P
        F_M = self.F_M
        baseline_purchase_costs = self.baseline_purchase_costs
        purchase_costs = self.purchase_costs
        installed_costs = self.installed_costs
        parallel = self.parallel
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        integer = int
        N_default = integer(parallel.get('self', 1))
        
        # Load main costs
        if N_default != 1:
            heat_utilities.extend((N_default - 1) * heat_utilities)
            power_utility.scale(N_default)
        for i in purchase_costs:
            if i not in baseline_purchase_costs:
                warning = RuntimeWarning(
                    f"Unit {self.ID}, adding items to the `purchase_costs` dictionary is "
                    "deprecated; add items to `baseline_purchase_costs` "
                    "dictionary instead"
                 )
                warn(warning)
                baseline_purchase_costs[i] = purchase_costs[i]
        for name, Cpb in baseline_purchase_costs.items(): 
            N = integer(parallel.get(name, N_default))
            if N == 1:
                if name in installed_costs and name in purchase_costs:
                    continue # Assume costs already added elsewhere using another method
            else:
                Cpb *= N
                baseline_purchase_costs[name] = Cpb
                if name in installed_costs and name in purchase_costs:
                    purchase_costs[name] *= N
                    installed_costs[name] *= N
                    continue # Assume costs already added elsewhere using another method
            F = F_D.get(name, 1.) * F_P.get(name, 1.) * F_M.get(name, 1.)
            try:
                installed_costs[name] = Cpb * (F_BM[name] + F - 1.)
            except KeyError:
                warn(f"the purchase cost item, '{name}', has "
                      "no defined bare-module factor in the "
                     f"'{type(self).__name__}.F_BM' dictionary; "
                      "bare-module factor now has a default value of 1",
                      RuntimeWarning)
                F_BM[name] = 1.
                installed_costs[name] = purchase_costs[name] = Cpb * F
            else:
                purchase_costs[name] = Cpb * F
        
        # Load auxiliary costs
        isa = isinstance
        for name, unit in self.get_auxiliary_units_with_names():
            if isa(unit, Unit):
                if not (unit._design or unit._cost): continue
            unit._load_costs() # Just in case user did not simulate or run summary.
            pname, *_ = name.split('[')
            N = integer(parallel.get(pname, N_default))
            if N == 1:
                heat_utilities.extend(unit.heat_utilities)
                power_utility.consumption += unit.power_utility.consumption
                power_utility.production += unit.power_utility.production
            elif N > self.max_parallel_units:
                raise RuntimeError(f'cannot have over a {self.max_parallel_units} unit operations in parallel')
            else:
                heat_utilities.extend(N * unit.heat_utilities)
                power_utility.consumption += N * unit.power_utility.consumption
                power_utility.production += N * unit.power_utility.production
            bpc_auxiliary = unit.baseline_purchase_costs
            pc_auxiliary = unit.purchase_costs
            ic_auxiliary = unit.installed_costs
            for i in bpc_auxiliary:
                j = ' - '.join([name.capitalize().replace('_', ' '), i])
                if j in baseline_purchase_costs: 
                    raise RuntimeError(
                        f"'{j}' already in `baseline_purchase_cost` "
                        f"dictionary of {repr(self)}; try using a different key"
                    )
                elif N == 1:
                    baseline_purchase_costs[j] = bpc_auxiliary[i]
                    purchase_costs[j] = pc_auxiliary[i]
                    installed_costs[j] = ic_auxiliary[i]
                else:
                    baseline_purchase_costs[j] = N * bpc_auxiliary[i]
                    purchase_costs[j] = N * pc_auxiliary[i]
                    installed_costs[j] = N * ic_auxiliary[i]
            
        self._costs_loaded = True
    
    def _setup(self):
        """Clear all results, setup up stream conditions and constant data, 
        and update system configuration based on units impacted by process 
        specifications. This method is run at the start of unit simulation, 
        before running mass and energy balances."""
        self.materialize_connections()
        self.power_utility.empty()
        for i in self.heat_utilities: i.empty()
        if not hasattr(self, '_N_heat_utilities'): self.heat_utilities.clear()
        self.parallel.clear()
        self.baseline_purchase_costs.clear()
        self.purchase_costs.clear()
        self.installed_costs.clear()
        self._costs_loaded = False
        for i in self.auxiliary_units: i._setup()
    
    def _check_setup(self):
        if any([self.power_utility, 
                self.heat_utilities, 
                self.baseline_purchase_costs, 
                self.purchase_costs, 
                self.installed_costs]):
            raise UnitInheritanceError(
               f'`{type(self).__name__}._setup` method did not clear unit results; '
                'a potential solution is to run `super()._setup()` in the `_setup` '
                'method of the unit subclass'
            )
    
    def _check_run(self):
        if any([self.power_utility, 
                self.heat_utilities, 
                self.baseline_purchase_costs, 
                self.purchase_costs, 
                self.installed_costs]):
            warn(
               f'`{type(self).__name__}._run` method added unit results '
                '(e.g., purchase costs, heat and power utilities); unit results '
                'should only be added in `_design` or `_cost` methods',
                RuntimeWarning
            )
    
    def materialize_connections(self):
        for s in self._ins + self._outs: 
            if not s: s.materialize_connection()
    
    @property
    def owner(self) -> Unit:
        owner = getattr(self, '_owner', None)
        return self if owner is None else owner.owner
    @owner.setter
    def owner(self, owner):
        if owner is not self: self._owner = owner
    
    def _get_tooltip_string(self, format=None, full=None):
        """Return a string that can be used as a Tippy tooltip in HTML output"""
        if format is None: format = bst.preferences.graphviz_format
        if full is None: full = bst.preferences.tooltips_full_results
        if format not in ('html', 'svg'): return ''
        if format == 'html' and full:
            results = self.results(include_installed_cost=True)
            tooltip = (
                " " + # makes sure graphviz does not try to parse the string as HTML
                results.to_html(justify='unset'). # unset makes sure that table header style can be overwritten in CSS
                replace("\n", "").replace("  ", "") # makes sure tippy.js does not add any whitespaces
            )
        else:
            newline = '<br>' if format == 'html' else '\n'
            electricity_consumption = self.power_utility.consumption
            electricity_production = self.power_utility.production
            cooling = self.net_cooling_duty / 1e3
            heating = self.net_heating_duty / 1e3
            utility_cost = self.utility_cost
            purchase_cost = int(float(self.purchase_cost))
            installed_cost = int(float(self.installed_cost))
            tooltip = ''
            if electricity_consumption:
                tooltip += f"{newline}Electricity consumption: {electricity_consumption:.3g} kW"
            if electricity_production:
                tooltip += f"{newline}Electricity production: {electricity_production:.3g} kW"
            if cooling:
                tooltip += f"{newline}Cooling duty: {cooling:.3g} MJ/hr"
            if heating:
                tooltip += f"{newline}Heating duty: {heating:.3g} MJ/hr"
            if utility_cost:
                tooltip += f"{newline}Utility cost: {utility_cost:.3g} USD/hr"
            if purchase_cost:
                tooltip += f"{newline}Purchase cost: {purchase_cost:,} USD"
            if installed_cost:
                tooltip += f"{newline}Installed equipment cost: {installed_cost:,} USD"
            if not tooltip: tooltip = 'No capital costs or utilities'
            elif tooltip: tooltip = tooltip.lstrip(newline)
            if format == 'html': tooltip = ' ' + tooltip
        return tooltip
    
    def get_design_input_units(self, key):
        design_input = getattr(self.__class__, key)
        return design_input.units
    
    def get_design_input(self, key: str, units: str):
        """
        Return design input in a new set of units of measure.
            
        Parameters
        ----------
        key :
            Name of design input.
        units :
            Units of measure.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed, tau=672)
        >>> tank.get_design_input('tau', 'day')
        28.0
        
        """
        return convert(getattr(self, key), self.get_design_input_units(key), units)
    
    def set_design_input(self, key: str, units: str, value: float):
        """
        Set design result in given the units of measure.
            
        Parameters
        ----------
        key :
            Name of design input.
        units :
            Units of measure.
        value:
            Value of the design input.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed)
        >>> tank.set_design_input('tau', 'day', 7)
        168
        >>> tank.get_design_input('tau', 'hr')
        168
        
        """
        value = convert(value, units, self.get_design_input_units(key))
        setattr(self, key, value)
        return value
    
    def convert_design_input(self, key, units, value):
        """
        Convert design input in given units to the stored units of measure.
            
        Parameters
        ----------
        key :
            Name of design input.
        units :
            Units of measure.
        value:
            Value of the design input.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed)
        >>> tank.simulate()
        >>> tank.convert_design_input('tau', 'day', 7)
        168
        
        """
        return convert(value, units, self.get_design_input_units(key))
    
    def get_design_result(self, key: str, units: str):
        """
        Return design result in a new set of units of measure.
            
        Parameters
        ----------
        key :
            Name of design result.
        units :
            Units of measure.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed)
        >>> tank.simulate()
        >>> tank.get_design_result('Total volume', 'm3')
        1214.19
        >>> tank.get_design_result('Total volume', 'L')
        1214191.0
        
        """
        return convert(self.design_results[key], self._units[key], units)
    
    def set_design_result(self, key: str, units: str, value: float):
        """
        Set design result in given the units of measure.
            
        Parameters
        ----------
        key :
            Name of design result.
        units :
            Units of measure.
        value:
            Value of the design result.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed)
        >>> tank.simulate()
        >>> tank.set_design_result('Total volume', 'm3', 1000)
        1000
        >>> tank.get_design_result('Total volume', 'm3')
        1000.0
        
        """
        self.design_results[key] = value = convert(value, units, self._units[key])
        return value
    
    def convert_design_result(self, key, units, value):
        """
        Convert design result in given units to the stored units of measure.
            
        Parameters
        ----------
        key :
            Name of design result.
        units :
            Units of measure.
        value:
            Value of the design result.
        
        Examples
        --------
        >>> import biosteam as bst
        >>> bst.settings.set_thermo(['Water'], cache=True)
        >>> feed = bst.Stream(None, Water=100)
        >>> tank = bst.StorageTank(None, feed)
        >>> tank.simulate()
        >>> tank.convert_design_result('Total volume', 'ft3', 1000)
        28.31
        
        """
        return convert(value, units, self._units[key])
    
    def reset_cache(self, isdynamic=None):
        pass
    
    def _get_design_input_info(self): 
        cls = self.__class__
        design_inputs = []
        for key in self._design_inputs:
            value = getattr(self, key)
            if value is None: continue # Decision not made
            dd = getattr(cls, key)
            design_inputs.append(
                (dd.name, value, dd.units)
            )
        return design_inputs
    
    def _load_stream_links(self):
        if self._link_streams: self._outs[0].link_with(self._ins[0])
    
    def _reevaluate(self):
        """Reevaluate design/cost/LCA results."""
        self._setup()
        self._summary()
    
    def _check_utilities(self):
        auxiliary_heat_utilities = set(sum([i.heat_utilities for i in self.auxiliary_units], []))
        for i in self.heat_utilities:
            if i in auxiliary_heat_utilities:
                raise UnitInheritanceError(
                    'auxiliary heat utilities were manually added to main utilities; '
                    'note that utilities from auxiliary units are already automatically '
                    'added to main unit operation'
                )
    
    def _summary(self, design_kwargs=None, cost_kwargs=None, lca_kwargs=None):
        """Run design/cost/LCA algorithms and compile results."""
        self._check_run()
        if not (self._design or self._cost): return
        if not self._skip_simulation_when_inlets_are_empty or not all([i.isempty() for i in self._ins]): 
            self._design(**design_kwargs) if design_kwargs else self._design()
            self._cost(**cost_kwargs) if cost_kwargs else self._cost()
            self._lca(**lca_kwargs) if lca_kwargs else self._lca()
            self._check_utilities()
        self._load_costs()
        self._load_operation_costs()

    def _load_operation_costs(self):
        ins = self._ins._streams
        outs = self._outs._streams
        prices = bst.stream_prices
        self._utility_cost = (
            sum([i.cost for i in self.heat_utilities]) 
            + self.power_utility.cost
            + sum([s.F_mass * prices[name] for name, index in self._inlet_utility_indices.items() if (s:=ins[index]).price == 0.])
            - sum([s.F_mass * prices[name] for name, index in self._outlet_utility_indices.items() if (s:=outs[index]).price == 0.])
        )
        fractions = self._get_flow_fractions()
        self._inlet_cost = sum(
            [(fractions[name] * ins[index].F_mass * prices[name]
              if name in fractions 
              else ins[index].F_mass * prices[name])
             for name, index in self._inlet_cost_indices.items()]
        )
        self._outlet_revenue = sum(
            [(fractions[name] * outs[index].F_mass * prices[name]
              if name in fractions 
              else outs[index].F_mass * prices[name])
             for name, index in self._outlet_revenue_indices.items()]
        )
    
    @property
    def specifications(self) -> list[ProcessSpecification]:
        """
        Process specifications as a list of process specification objects.
        
        See Also
        --------
        add_specification
        add_bounded_numerical_specification
        
        """
        return self._specifications
    @specifications.setter
    def specifications(self, specifications):
        if specifications is None:
            self._specifications = []
        else:
            self._specifications = specifications
    
    @property
    def specification(self):
        raise AttributeError('`specification` property is deprecated; use `add_specification` or `specifications` (plural with an s) instead')
    @specification.setter
    def specification(self, specification):
        raise AttributeError('`specification` property is deprecated; use `add_specification` or `specifications` (plural with an s) instead')
    
    @property
    def baseline_purchase_cost(self) -> float:
        """Total baseline purchase cost, without accounting for design ,
        pressure, and material factors [USD]."""
        return sum(self.baseline_purchase_costs.values())
    
    @property
    def purchase_cost(self) -> float:
        """Total purchase cost [USD]."""
        return sum(self.purchase_costs.values())
    
    @property
    def installed_cost(self) -> float:
        """Total installed equipment cost [USD]."""
        return sum(self.installed_costs.values())
    
    @property
    def utility_cost(self) -> float:
        """Total utility cost [USD/hr]."""
        return self._utility_cost

    def mass_balance_error(self):
        """Return error in stoichiometric mass balance. If positive,
        mass is being created. If negative, mass is being destroyed."""
        return self.F_mass_out - self.F_mass_in

    def atomic_balance_error(self):
        """Return a dictionary of errors in stoichiometric atomic balances. 
        If value is positive, the atom is being created. If negative, the atom 
        is being destroyed."""
        from chemicals import elements
        mol = sum([i.mol for i in self._outs._streams]) - sum([i.mol for i in self._ins._streams])
        formula_array = self.chemicals.formula_array
        unbalanced_array = formula_array @ mol
        return elements.array_to_atoms(unbalanced_array)

    def empty(self):
        """
        Empty all unit operation results and outlet flows.
        """
        self.design_results.clear()
        self.baseline_purchase_costs.clear()
        self.purchase_costs.clear()
        self.installed_costs.clear()
        for i in self.outs: i.empty()
        self.heat_utilities.clear()
        self.power_utility.empty()
        self._utility_cost = 0.
        
    def simulate(self, 
            run: Optional[bool]=None,
            design_kwargs: Optional[dict]=None,
            cost_kwargs: Optional[dict]=None):
        """
        Run rigorous simulation and determine all design requirements.
        
        Parameters
        ----------
        run :
            Whether to run mass and energy balance or to assume the same inlet
            and outlet conditions. Defaults to True.
        design_kwargs :
            Keyword arguments passed to `_design` method.
        cost_kwargs :
            Keyword arguments passed to `_cost` method.
            
        """
        self._setup()
        self._check_setup()
        if run is None or run:
            for ps in self._specifications: ps.compile_path(self)
            self._load_stream_links()
            self.run()
        self._summary(design_kwargs, cost_kwargs)

    def results(self, 
            with_units=True, 
            include_utilities=True,
            include_total_cost=True, 
            include_installed_cost=True,
            include_zeros=True, 
            external_utilities=None, 
            key_hook=None,
            basis=None
        ):
        """
        Return key results from simulation as a DataFrame if `with_units`
        is True or as a Series otherwise.
        """
        def addkey(key):
            if key_hook: key = key_hook(key)
            keys.append(key)
            
        def addcapex(key):
            if key_hook: key = key_hook(key)
            *others, name = key
            names = name.split(' - ')
            parent = self
            parallel = parent.parallel
            N = N_default = parallel.get('self', 1)
            for auxname in names:
                auxsearch = auxname.replace(' ', '_').lower()
                N *= int(parallel.get(auxname, 1.))
                for i in parent.auxiliary_unit_names:
                    if auxsearch == i.lower():
                        parent = getattr(parent, i)
                        if isinstance(parent, list):
                            N *= len(parent)
                        else:
                            parallel = parent.parallel
                            N_default = parallel.get('self', 1)
                            N *= int(parallel.get(i, N_default))
                        break
                else:
                    break
            if N != 1: key = (*others, name + f' (x{N})')
            keys.append(key)
        keys = []; 
        vals = []; addval = vals.append
        stream_prices = bst.stream_prices
        all_utilities = self.heat_utilities + external_utilities if external_utilities else self.heat_utilities
        if with_units:
            if include_utilities:
                power_utility = self.power_utility
                if power_utility:
                    addkey(('Electricity', 'Power'))
                    addval(('kW', power_utility.power))
                    if include_zeros or power_utility.cost:
                        addkey(('Electricity', 'Cost'))
                        addval(('USD/hr', power_utility.cost))
                for heat_utility in HeatUtility.sum_by_agent(all_utilities):
                    if heat_utility:
                        ID = heat_utility.ID.replace('_', ' ').capitalize()
                        addkey((ID, 'Duty'))
                        addval(('kJ/hr', heat_utility.duty))
                        addkey((ID, 'Flow'))
                        addval(('kmol/hr', heat_utility.flow))
                        if include_zeros or heat_utility.cost: 
                            addkey((ID, 'Cost'))
                            addval(('USD/hr', heat_utility.cost))
                for name, flow in self.get_inlet_cost_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (inlet)'
                        addkey((ID, 'Flow'))
                        addval(('kg/hr', flow))
                        addkey((ID, 'Cost'))
                        addval(('USD/hr', flow * stream_prices[name]))
                for name, flow in self.get_outlet_revenue_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (outlet)'
                        addkey((ID, 'Flow'))
                        addval(('kg/hr', flow))
                        addkey((ID, 'Cost'))
                        addval(('USD/hr', - flow * stream_prices[name]))
            units = self._units
            Cost = self.purchase_costs
            if basis:
                ureg = tmo.units_of_measure.ureg
                ureg.default_system = basis
                Quantity = ureg.Quantity
                def to_basis(vi, ui):
                    if not ui or ui == 'hr' or isinstance(vi, str): return ui, vi
                    try:
                        q = Quantity(vi, ui).to_base_units()
                    except Exception: # Units not defined in registry (e.g., psig)
                        return ui, vi
                    return q.units, q.magnitude
                for ki, vi, ui in self._get_design_input_info():
                    addkey(('Design inputs', ki))
                    addval(to_basis(vi, ui))
                for ki, vi in self.design_results.items():
                    addkey(('Design results', ki))
                    addval(to_basis(vi, units.get(ki, '')))
            else:
                for ki, vi, ui in self._get_design_input_info():
                    addkey(('Design inputs', ki))
                    addval((ui, vi))
                for ki, vi in self.design_results.items():
                    addkey(('Design results', ki))
                    addval((units.get(ki, ''), vi))
            for ki, vi in Cost.items():
                addcapex(('Purchase cost', ki))
                addval(('USD', vi))
            if include_total_cost:
                addkey(('Total purchase cost', ''))
                addval(('USD', self.purchase_cost))
                if include_installed_cost:
                    addkey(('Installed equipment cost', ''))
                    addval(('USD', self.installed_cost))
                utility_cost = self.utility_cost
                if include_zeros or utility_cost: 
                    addkey(('Utility cost', ''))
                    addval(('USD/hr', utility_cost))
            if not keys: return None
            df = pd.DataFrame(vals,
                              pd.MultiIndex.from_tuples(keys),
                              ('Units', self.ID))
            df.columns.name = self.line
            return df
        else:
            if include_utilities:
                power_utility = self.power_utility
                if power_utility:
                    addkey(('Electricity', 'Power'))
                    addval(power_utility.power)
                    if include_zeros or power_utility.cost:
                        addkey(('Electricity', 'Cost'))
                        addval(power_utility.cost)
                for heat_utility in HeatUtility.sum_by_agent(all_utilities):
                    if heat_utility:
                        ID = heat_utility.ID.replace('_', ' ').capitalize()
                        addkey((ID, 'Duty'))
                        addval(heat_utility.duty)
                        addkey((ID, 'Flow'))
                        addval(heat_utility.flow)
                        if include_zeros or heat_utility.cost:
                            addkey((ID, 'Cost'))
                            addval(heat_utility.cost)
                for name, flow in self.get_inlet_cost_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (inlet)'
                        addkey((ID, 'Flow'))
                        addval(flow)
                        addkey((ID, 'Cost'))
                        addval(flow * stream_prices[name])
                for name, flow in self.get_outlet_revenue_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (outlet)'
                        addkey((ID, 'Flow'))
                        addval(flow)
                        addkey((ID, 'Cost'))
                        addval(-flow * stream_prices[name])
            for ki, vi, ui in self._get_design_input_info():
                addkey(('Design inputs', ki))
                addval(vi)
            for ki, vi in self.design_results.items():
                addkey(('Design results', ki))
                addval(vi)
            for ki, vi in self.purchase_costs.items():
                addcapex(('Purchase cost', ki))
                addval(vi)
            if include_total_cost:
                addkey(('Total purchase cost', ''))
                addval(self.purchase_cost)
                if include_installed_cost:
                    addkey(('Installed equipment cost', ''))
                    addval(self.installed_cost)
                utility_cost = self.utility_cost
                if include_zeros or utility_cost:
                    addkey(('Utility cost', ''))
                    addval(utility_cost)
            if not keys: return None
            series = pd.Series(vals, pd.MultiIndex.from_tuples(keys))
            series.name = self.ID
            return series

    def get_available_chemicals(self):
        streams = [i for i in (self._ins + self._outs) if i]
        reaction_chemicals = sum([i.reaction_chemicals for i in self.__dict__.values() if hasattr(i, 'reaction_chemicals')], [])
        required_chemicals = set(sum([i.available_chemicals for i in streams], reaction_chemicals))
        return [i for i in self.chemicals if i in required_chemicals]

    ### Net input and output flows ###
    
    # Molar flow rates
    @property
    def mol_in(self) -> NDArray[float]:
        """Molar flows going in [kmol/hr]."""
        return sum([s.mol for s in self._ins if s])
    @property
    def mol_out(self) -> NDArray[float]:
        """Molar flows going out [kmol/hr]."""
        return sum([s.mol for s in self._outs if s])

    @property
    def z_mol_in(self) -> NDArray[float]:
        """Molar fractions going in [kmol/hr]."""
        return self.mol_in/self.F_mol_in
    @property
    def z_mol_out(self) -> NDArray[float]:
        """Molar fractions going in."""
        return self.mol_out/self.F_mol_out

    @property
    def F_mol_in(self) -> float:
        """Net molar flow going in [kmol/hr]."""
        return sum([s.F_mol for s in self._ins if s])
    @property
    def F_mol_out(self) -> float:
        """Net molar flow going out [kmol/hr]."""
        return sum([s.F_mol for s in self._outs if s])

    # Mass flow rates
    @property
    def mass_in(self)-> NDArray[float]:
        """Mass flows going in [kg/hr]."""
        return sum([s.mol for s in self._ins if s]) * self._thermo.chemicals.MW
    @property
    def mass_out(self)-> NDArray[float]:
        """Mass flows going out [kg/hr]."""
        return sum([s.mol for s in self._outs if s]) * self._thermo.chemicals.MW

    @property
    def z_mass_in(self)-> NDArray[float]:
        """Mass fractions going in."""
        return self.mass_in/self.F_mass_in
    @property
    def z_mass_out(self)-> NDArray[float]:
        """Mass fractions going out."""
        return self.mass_out/self.F_mass_out

    @property
    def F_mass_in(self)-> float:
        """Net mass flow going in [kg/hr]."""
        return self.mass_in.sum()
    @property
    def F_mass_out(self) -> float:
        """Net mass flow going out [kg/hr]."""
        return self.mass_out.sum()

    # Volumetric flow rates
    @property
    def vol_in(self) -> NDArray[float]:
        """Volumetric flows going in [m3/hr]."""
        return sum([s.vol for s in self._ins if s])
    @property
    def F_vol_in(self) -> float:
        """Net volumetric flow going in [m3/hr]."""
        return sum([s.F_vol for s in self._ins if s])

    @property
    def z_vol_in(self) -> NDArray[float]:
        """Volumetric fractions going in."""
        return self.vol_in/self.F_vol_in
    @property
    def vol_out(self) -> NDArray[float]:
        """Volumetric flows going out [m3/hr]."""
        return sum([s.vol for s in self._outs if s])

    @property
    def F_vol_out(self)-> float:
        """Net volumetric flow going out [m3/hr]."""
        return sum([s.F_vol for s in self._outs if s])
    @property
    def z_vol_out(self) -> NDArray[float]:
        """Volumetric fractions going out."""
        return self.vol_out/self.F_vol_out

    # Enthalpy flow rates
    @property
    def H_in(self) -> float:
        """Enthalpy flow going in [kJ/hr]."""
        return sum([s.H for s in self._ins if s])

    @property
    def H_out(self) -> float:
        """Enthalpy flow going out [kJ/hr]."""
        return sum([s.H for s in self._outs if s])

    @property
    def Hf_in(self) -> float:
        """Enthalpy of formation flow going in [kJ/hr]."""
        return sum([s.Hf for s in self._ins if s])

    @property
    def Hf_out(self) -> float:
        """Enthalpy of formation flow going out [kJ/hr]."""
        return sum([s.Hf for s in self._outs if s])

    @property
    def Hnet(self) -> float:
        """Net enthalpy flow, including enthalpies of formation [kJ/hr]."""
        return self.H_out - self.H_in + self.Hf_out - self.Hf_in


# %% Process specification

class ProcessSpecification:
    __slots__ = ('f', 'args', 'impacted_units', 'path')
    
    def __init__(self, f, args, impacted_units):
        self.f = f
        self.args = args
        if impacted_units: 
            self.impacted_units = tuple(impacted_units)
        else:
            self.impacted_units = self.path = ()
        
    def __call__(self):
        self.f(*self.args)
        isa = isinstance
        for i in self.path: 
            if isa(i, Unit): i.run()
            else: i.converge() # Must be a system
            
    def compile_path(self, unit):
        self.path = unit.path_from(self.impacted_units, system=unit._system) if self.impacted_units else ()
        
    def create_temporary_connections(self, unit):
        # Temporary connections are created first than the path because 
        # temporary connections may change the system configuration 
        # such that the path is incorrect.
        impacted_units = self.impacted_units
        if impacted_units:
            upstream_units = unit.get_upstream_units()
            for other in impacted_units:
                if other not in upstream_units: temporary_connection(unit, other)
            
    def __repr__(self):
        return f"{type(self).__name__}(f={display_asfunctor(self.f)}, args={self.args}, impacted_units={self.impacted_units})"


class BoundedNumericalSpecification:
    __slots__ = (
        'f', 'x0', 'x1', 'y0', 'y1', 'x', 'xtol', 'ytol', 'args', 
        'maxiter', 'checkroot', 'checkiter', 'checkbounds', 'x_last',
    )
    
    def __init__(self, 
            f: Callable,
            x0: float, 
            x1: float, 
            y0: Optional[float]=None, 
            y1: Optional[float]=None, 
            x: Optional[float]=None,
            xtol: float=0.,
            ytol: float=5e-8,
            args: Tuple[Any, ...]=(), 
            maxiter: int=50,
            checkroot: bool=False, 
            checkiter: bool=True, 
            checkbounds: bool=True
        ):
        self.f = f
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.x = x
        self.xtol = xtol
        self.ytol = ytol
        self.args = args
        self.maxiter = maxiter
        self.checkroot = checkroot
        self.checkiter = checkiter
        self.checkbounds = checkbounds
        self.x_last = None
        
    def __call__(self):
        self.x = IQ_interpolation(
            self.f, self.x0, self.x1, self.y0, self.y1, self.x, self.xtol, self.ytol, 
            self.args, self.maxiter, self.checkroot, self.checkiter, self.checkbounds,
        )
        return self.x
    
    def compile_path(self, unit): pass
    
    def create_temporary_connections(self, unit): pass


# %% Temporary units for temporary connections between unit operations
#    (so that process specifications with impacted units are simulated in the right order)

class TemporaryStream(Stream):
    """Empty stream used only to create temporary connections between units."""
    __slots__ = ()
    
    def __init__(self, ID=None, source=None, sink=None, thermo=None):
        self._register(ID)
        self._source = source
        self._sink = sink
        self._load_thermo(thermo)


temporary_units_dump = []

def temporary_connection(source, sink):
    upstream = source.outs[0]
    downstream = sink.ins[0]
    thermo = source.thermo
    temporary_stream = TemporaryStream(None, thermo=thermo)
    if isinstance(upstream.sink, TemporarySource):
        upstream.sink.outs.append(temporary_stream)
    else:
        stream = TemporaryStream(None, thermo=thermo)
        old_connection = upstream.get_connection()
        sink = upstream.sink
        if sink: upstream.sink.ins.replace(upstream, stream)
        TemporarySource(upstream, [stream, temporary_stream], old_connection)
    if isinstance(downstream.source, TemporarySink):
        downstream.source.ins.append(temporary_stream)
    else:
        stream = TemporaryStream(None, thermo=thermo)
        old_connection = downstream.get_connection()
        downstream.sink.ins.replace(downstream, stream)
        TemporarySink([downstream, temporary_stream], stream, old_connection)


class TemporaryUnit:
    __slots__ = ('ID', '_ID', 'ins', 'outs', '_ins', '_outs', 'old_connection')
    _universal = False
    _interaction = False
    auxiliary_units = ()
    def __init__(self, ins, outs, old_connection):
        temporary_units_dump.append(self)
        self.ID = self._ID = 'TU' + str(len(temporary_units_dump))
        self.old_connection = old_connection
        self.ins = self._ins = piping.Inlets(
            self, self._N_ins, ins, None, self._ins_size_is_fixed, 5
        )
        self.outs = self._outs = piping.Outlets(
            self, self._N_outs, outs, None, self._outs_size_is_fixed, 5
        )
        
    neighborhood = Unit.neighborhood
    get_downstream_units = Unit.get_downstream_units
    get_upstream_units = Unit.get_upstream_units
    _add_upstream_neighbors_to_set = Unit._add_upstream_neighbors_to_set
    _add_downstream_neighbors_to_set = Unit._add_downstream_neighbors_to_set
    __repr__ = Unit.__repr__
    __str__ = Unit.__str__
        

class TemporarySource(TemporaryUnit):
    __slots__ = ()
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = True
    _outs_size_is_fixed = False
    

class TemporarySink(TemporaryUnit):
    __slots__ = ()
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True


# %% Design and capital

class UnitDesignAndCapital:
    
    __slots__ = (
        'unit',
        'F_BM',
        'F_D',
        'F_P',
        'F_M',
        'design_results',
        'baseline_purchase_costs',
        'purchase_costs',
        'installed_costs',
    )
    
    baseline_purchase_cost = Unit.baseline_purchase_cost
    purchase_cost = Unit.purchase_cost
    installed_cost = Unit.installed_cost
    
    def __init__(self, 
            unit, F_BM, F_D, F_P, F_M,
            design_results: dict,
            baseline_purchase_costs: dict,
            purchase_costs: dict,
            installed_costs: dict,
        ):
        self.unit = unit
        self.F_BM = F_BM
        self.F_D = F_D
        self.F_P = F_P
        self.F_M = F_M
        self.design_results = design_results
        self.baseline_purchase_costs = baseline_purchase_costs
        self.purchase_costs = purchase_costs
        self.installed_costs = installed_costs

    @property
    def equipment_lifetime(self):
        return self.unit.equipment_lifetime
