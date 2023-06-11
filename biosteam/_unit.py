# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from __future__ import annotations
from .utils import ignore_docking_warnings
import numpy as np
import pandas as pd
from warnings import warn
from ._graphics import UnitGraphics, box_graphics
from ._heat_utility import HeatUtility
from .utils import AbstractMethod, format_title, static, piping, StreamLinkOptions
from ._power_utility import PowerUtility
from .exceptions import UnitInheritanceError
from thermosteam.utils import thermo_user, registered
from thermosteam.units_of_measure import convert
from copy import copy
import biosteam as bst
from thermosteam import Stream
from thermosteam.base import display_asfunctor
from numpy.typing import NDArray
from typing import Callable, Optional, TYPE_CHECKING, Sequence, Iterable
import thermosteam as tmo
if TYPE_CHECKING: 
    System = bst.System
    HXutility = bst.HXutility
    UtilityAgent = bst.UtilityAgent

streams = Optional[Sequence[Stream]] # TODO: Replace with Stream|str and make it an explicit TypeAlias once BioSTEAM moves to Python 3.10
int_types = (int, np.int32)

__all__ = ('Unit',)

_count = [0]
def count():
    _count[0] += 1
    print(_count)

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
        impacted_units = self.impacted_units
        self.path = unit.path_from(impacted_units, system=unit._system) if impacted_units else ()
        
    def create_temporary_connections(self, unit):
        # Temporary connections are created first than the path because 
        # temporary connections may change the system configuration 
        # such that the path is incorrect.
        impacted_units = self.impacted_units
        if impacted_units:
            downstream_units = unit.get_downstream_units()
            upstream_units = unit.get_upstream_units()
            connected_units = upstream_units | downstream_units
            for other in impacted_units:
                if other not in connected_units:
                    bst.temporary_connection(unit, other)
            
    def __repr__(self):
        return f"{type(self).__name__}(f={display_asfunctor(self.f)}, args={self.args}, impacted_units={self.impacted_units})"

# %% Typing

# from typing import Collection, Union, Annotated
# streams = Union[Collection[Union[Stream, str, None]], Union[Stream, str, None]]
# stream = Union[Annotated[Union[Stream, str, None], 1], Union[Stream, str, None]]
# stream_sequence = Collection[Union[Stream, str, None]]

# %% Inlet and outlet representation

def repr_ins_and_outs(layout, ins, outs, T, P, flow, composition, N, IDs, sort, data):
    info = ''
    if ins:
        info += 'ins...\n'
        i = 0
        for stream in ins:
            unit = stream._source
            source_info = f'  from  {type(unit).__name__}-{unit}' if unit else ''
            if stream and data:
                stream_info = stream._info(layout, T, P, flow, composition, N, IDs, sort)
                index = stream_info.index('\n')
                number = f'[{i}] '
                spaces = len(number) * ' '
                info += number + str(stream) + source_info + stream_info[index:].replace('\n', '\n' + spaces) + '\n'
            else:
                info += f'[{i}] {stream}' + source_info + '\n'
            i += 1
    if outs:
        info += 'outs...\n'
        i = 0
        for stream in outs:
            unit = stream._sink
            sink_info = f'  to  {type(unit).__name__}-{unit}' if unit else ''
            if stream and data:
                stream_info = stream._info(layout, T, P, flow, composition, N, IDs, sort)
                index = stream_info.index('\n')
                number = f'[{i}] '
                spaces = len(number) * ' '
                info += number + str(stream) + sink_info + stream_info[index:].replace('\n', '\n' + spaces) + '\n'
            else:
                info += f'[{i}] {stream}' + sink_info + '\n'
            i += 1
    return info[:-1]

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
    
# %% Unit Operation

@thermo_user
@registered(ticket_name='U')
class Unit:
    """
    Abstract class for Unit objects. Child objects must contain
    :attr:`~Unit._run`, :attr:`~Unit._design` and :attr:`~Unit._cost` methods to 
    estimate stream outlets of a Unit and find design and cost information.  

    Parameters
    ----------
    ID :
        A unique identification. If ID is None, unit will not be
        registered in flowsheet. By default, a unique ID will be chosen.
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
    
    def __init_subclass__(cls,
                          isabstract=False,
                          new_graphics=True,
                          does_nothing=None):
        if does_nothing: return 
        dct = cls.__dict__
        if '_N_heat_utilities' in dct:
            warn("'_N_heat_utilities' class attribute is scheduled for deprecation; "
                 "use the `add_heat_utility` method instead",
                 DeprecationWarning, stacklevel=2)
        if 'run' in dct:
            raise UnitInheritanceError(
                 "the 'run' method cannot be overridden; implement `_run` instead"
            )
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
            if hasattr(cls, '_BM'): 
                raise UnitInheritanceError(
                    'cannot set `_BM`; implement `_F_BM_default` instead'
                )
            elif hasattr(cls, '_F_BM_defaults'):
                raise UnitInheritanceError(
                    'cannot set `_F_BM_defaults`; implement '
                    '`_F_BM_default` instead'
                )
            elif cls._F_BM_default is Unit._F_BM_default:
                cls._F_BM_default = {}
            
            if hasattr(cls, '_equipment_lifetime'):
                raise UnitInheritanceError(
                    'cannot set `_equipment_lifetime`; '
                    'implement `_default_equipment_lifetime` instead'
                )
            elif hasattr(cls, '_default_equipment_lifetimes'):
                raise UnitInheritanceError(
                    'cannot set `_default_equipment_lifetimes`; implement '
                    '`_default_equipment_lifetime` instead'
                )
            elif cls._default_equipment_lifetime is Unit._default_equipment_lifetime: 
                cls._default_equipment_lifetime = {}
            if cls._units is Unit._units: cls._units = {}
            if not cls._run:
                if cls._N_ins == 1 and cls._N_outs == 1:
                    static(cls)
                else:
                    raise UnitInheritanceError(
                        "Unit subclass with multiple inlet or outlet streams "
                        "must implement a '_run' method unless the "
                        "'isabstract' keyword argument is True"
                    )
        if '__init__' in dct:
            init = dct['__init__']
            annotations = init.__annotations__
            for i in ('ins', 'outs'):
                if i not in annotations: annotations[i] = streams
            if '_stacklevel' not in dct: cls._stacklevel += 1
        name = cls.__name__
        if hasattr(bst, 'units') and hasattr(bst, 'wastewater') and hasattr(bst, 'facilities'):
            # Add 3rd party unit to biosteam module for convinience
            if name not in bst.units.__dict__:
                bst.units.__dict__[name] = cls
            if name not in bst.__dict__:
                bst.__dict__[name] = cls
        
    ### Abstract Attributes ###
    #: **class-attribute** Units of measure for :attr:`~Unit.design_results` dictionary.
    _units: dict[str, str] = {}
    
    #: **class-attribute** Default bare-module factors for each purchase cost item.
    #: Items in this dictionary are copied to the :attr:`~Unit.F_BM` attribute during 
    #: initialization.
    _F_BM_default: dict[str, float] = {}
    
    #: **class-attribute** Cost items that need to be summed across operation modes for 
    #: flexible operation (e.g., filtration membranes).
    _materials_and_maintenance: frozenset[str] = frozenset()
    
    #: **class-attribute** Name of attributes that are auxiliary units. These units
    #: will be accounted for in the purchase and installed equipment costs
    #: without having to add these costs in the :attr:`~Unit.baseline_purchase_costs` dictionary.
    #: Heat and power utilities are also automatically accounted for.
    auxiliary_unit_names: tuple[str, ...] = ()
    
    #: **class-attribute** Expected number of inlet streams. Defaults to 1.
    _N_ins: int = 1  
    
    #: **class-attribute** Expected number of outlet streams. Defaults to 1
    _N_outs: int = 1
    
    #: **class-attribute** Whether the number of streams in :attr:`~Unit.ins` is fixed.
    _ins_size_is_fixed: bool = True
    
    #: **class-attribute** Whether the number of streams in :attr:`~Unit.outs` is fixed.
    _outs_size_is_fixed: bool = True
    
    #: **class-attribute** Options for linking streams
    _stream_link_options: StreamLinkOptions = None

    #: **class-attribute** Used for piping warnings.
    _stacklevel: int = 5
    
    #: **class-attribute** Name denoting the type of Unit class. Defaults to the class
    #: name of the first child class
    line: str = 'Unit'

    #: **class-attribute** Lifetime of equipment. Defaults to lifetime of
    #: production venture. Use an integer to specify the lifetime for all
    #: items in the unit purchase costs. Use a dictionary to specify the 
    #: lifetime of each purchase cost item.
    _default_equipment_lifetime: int|dict[str, int] = {}

    #: **class-attribute** Settings for diagram representation. Defaults to a 
    #: box with the same number of inlet and outlet edges as :attr:`~Unit._N_ins` 
    #: and :attr:`~Unit._N_outs`.
    _graphics: UnitGraphics = box_graphics

    #: **class-attribute** Whether to skip detailed simulation when inlet 
    #: streams are empty. If inlets are empty and this flag is True,
    #: detailed mass and energy balance, design, and costing algorithms are skipped
    #: and all outlet streams are emptied.
    _skip_simulation_when_inlets_are_empty = False

    ### Abstract methods ###
    
    #: Create auxiliary components.
    _load_components = AbstractMethod
    
    #: Run mass and energy balances and update outlet streams (without user-defined specifications).
    _run = AbstractMethod
    
    #: Add design requirements to the :attr:`~Unit.design_results` dictionary.
    _design = AbstractMethod
    
    #: Add itemized purchase costs to the :attr:`~Unit.baseline_purchase_costs` dictionary.
    _cost = AbstractMethod

    #: For embodied emissions (e.g., unit construction) in LCA
    _lca = AbstractMethod
    
    def __init__(self, ID: Optional[str]='', ins=None, outs=(), thermo: tmo.Thermo=None):
        self._system = None
        self._isdynamic = False
        self._register(ID)
        self._load_thermo(thermo)
    
        ### Initialize streams
        
        self._ins = piping.Inlets(
            self, self._N_ins, ins, self._thermo, self._ins_size_is_fixed, self._stacklevel
        )
        self._outs = piping.Outlets(
            self, self._N_outs, outs, self._thermo, self._outs_size_is_fixed, self._stacklevel
        )
    
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
        
        try:
            #: Lifetime of equipment. Defaults to values in the class attribute 
            #: :attr:`~Unit._default_equipment_lifetime`. Use an integer to specify the lifetime 
            #: for all items in the unit purchase costs. Use a dictionary to specify 
            #: the lifetime of each purchase cost item.
            self.equipment_lifetime: int|dict[str, int] = copy(self._default_equipment_lifetime)
        except AttributeError:
            self.equipment_lifetime = {}
    
        ### Initialize specification    
    
        #: All specification functions
        self._specifications: list[Callable] = []
        
        #: Whether to run mass and energy balance after calling
        #: specification functions
        self.run_after_specifications: bool = False 
        
        #: Whether to prioritize unit operation specification within recycle loop (if any).
        self.prioritize: bool = False
        
        #: Safety toggle to prevent infinite recursion
        self._active_specifications: set[ProcessSpecification] = set()
        
        #: Name-number pairs of baseline purchase costs and auxiliary unit 
        #: operations in parallel. Use 'self' to refer to the main unit. Capital 
        #: and heat and power utilities in parallel will become proportional to this 
        #: value.
        self.parallel: dict[str, int] = {}
        
        self._assert_compatible_property_package()
        
        self._utility_cost = None
    
    def _init_ins(self, ins):
        self._ins = piping.Inlets(
            self, self._N_ins, ins, self._thermo, self._ins_size_is_fixed, self._stacklevel
        )
    
    def _init_outs(self, outs):
        self._outs = piping.Outlets(
            self, self._N_outs, outs, self._thermo, self._outs_size_is_fixed, self._stacklevel
        )

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
    
    @property
    def feed(self) -> Stream:
        """Equivalent to :attr:`~Unit.ins`\[0] when the number of inlets is 1."""
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
        """Equivalent to :attr:`~Unit.outs`\[0] when the number of outlets is 1."""
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
            Enforced fraction of heat transfered from utility (due
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
            Name of utility, as defined in :meth:`settings.stream_utility_prices <thermosteam._settings.ProcessSettings.stream_utility_prices>`.
        stream :
            Inlet or outlet utility stream.
        
        """
        if name not in bst.stream_utility_prices:
            raise ValueError(f"price of '{name}' must be defined in settings.stream_utility_prices")
        if stream._sink is self:
            self._inlet_utility_indices[name] = self._ins._streams.index(stream)
        elif stream._source is self:
            self._outlet_utility_indices[name] = self._outs._streams.index(stream)
        else:
            raise ValueError(f"stream '{stream.ID}' must be connected to {repr(self)}")
            
    def get_inlet_utility_flows(self):
        ins = self._ins._streams
        return {name: ins[index].F_mass for name, index in self._inlet_utility_indices.items()}
    
    def get_outlet_utility_flows(self):
        outs = self._outs._streams
        return {name: outs[index].F_mass for name, index in self._outlet_utility_indices.items()}
    
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
        if self._costs_loaded: return
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
            unit.owner = self # In case units are created dynamically
            if isa(unit, Unit):
                if not (unit._design or unit._cost): continue
            unit._load_costs() # Just in case user did not simulate or run summary.
            N = integer(parallel.get(name, N_default))
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
            F_BM_auxiliary = unit.F_BM
            F_D_auxiliary = unit.F_D
            F_P_auxiliary = unit.F_P
            F_M_auxiliary = unit.F_M
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
                else:
                    F_D[j] = fd = F_D_auxiliary.get(i, 1.)
                    F_P[j] = fp = F_P_auxiliary.get(i, 1.)
                    F_M[j] = fm = F_M_auxiliary.get(i, 1.)
                    if N == 1:
                        baseline_purchase_costs[j] = Cpb = bpc_auxiliary[i]
                        purchase_costs[j] = pc_auxiliary[i]
                        installed_costs[j] = Cbm = ic_auxiliary[i]
                    else:
                        baseline_purchase_costs[j] = Cpb = N * bpc_auxiliary[i]
                        purchase_costs[j] = N * pc_auxiliary[i]
                        installed_costs[j] = Cbm = N * ic_auxiliary[i]
                    try:
                        F_BM[j] = F_BM_auxiliary[i]
                    except KeyError:
                        # Assume costs already added elsewhere using another method.
                        # Calculate BM as an estimate.
                        F_BM[j] = Cbm / Cpb + 1 - fd * fp * fm
            
        self._costs_loaded = True
    
    def _setup(self):
        """Clear all results, setup up stream conditions and constant data, 
        and update system configuration based on units impacted by process 
        specifications. This method is run at the start of unit simulation, 
        before running mass and energy balances."""
        for i in self.auxiliary_units: i._setup()
        self.power_utility.empty()
        for i in self.heat_utilities: i.empty()
        if not hasattr(self, '_N_heat_utilities'): self.heat_utilities.clear()
        self.parallel.clear()
        self.baseline_purchase_costs.clear()
        self.purchase_costs.clear()
        self.installed_costs.clear()
        self._costs_loaded = False
    
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
            raise UnitInheritanceError(
               f'`{type(self).__name__}._run` method added unit results '
                '(e.g., purchase costs, heat and power utilities); unit results '
                'should only be added in `_design` or `_cost` methods'
            )
    
    def materialize_connections(self):
        for s in self._ins + self._outs: 
            if not s: s.materialize_connection()
    
    @property
    def system(self) -> System|None:
        return self._system
    
    @property
    def owner(self) -> Unit:
        owner = getattr(self, '_owner', None)
        return self if owner is None else owner.owner
    @owner.setter
    def owner(self, owner):
        self._owner = None if owner is self else owner
    
    @ignore_docking_warnings
    def disconnect(self, discard=False, inlets=None, outlets=None, join_ends=False):
        ins = self._ins
        outs = self._outs
        if inlets is None: 
            inlets = [i for i in ins if i.source]
            ins[:] = ()
        else:
            for i in inlets: ins[ins.index(i) if isinstance(i, Stream) else i] = None
        if outlets is None: 
            outlets = [i for i in outs if i.sink]
            outs[:] = ()
        else:
           for o in outlets: outs[ins.index(o) if isinstance(o, Stream) else o] = None
        if join_ends:
            if len(inlets) != len(outlets):
                raise ValueError("number of inlets must match number of outlets to join ends")
            for inlet, outlet in zip(inlets, outlets):
                outlet.sink.ins.replace(outlet, inlet)
        if discard: bst.main_flowsheet.discard(self)

    @ignore_docking_warnings
    def insert(self, stream, inlet=None, outlet=None):
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
        added_stream = False
        if outlet is None:
            if self._outs_size_is_fixed:
                if self._N_outs == 1:
                    sink.ins.replace(stream, self.outs[0])
                else:
                    raise ValueError("undefined outlet; must pass outlet when outlets are fixed and multiple are available")
            else:
                self.outs.append(stream)
                added_stream = True
        else:
            if isinstance(outlet, Stream) and outlet.source is not self:
                raise ValueError("source of given outlet must be this unit")
            else:
                outlet = self.outs[outlet]
            sink.ins.replace(stream, outlet)
        if inlet is None:
            if self._ins_size_is_fixed or added_stream:
                if self._N_ins == 1:
                    source.outs.replace(stream, self.ins[0])
                else:
                    raise ValueError("undefined inlet; must pass inlet when inlets are fixed and multiple are available")
            else:
                self.ins.append(stream)
        else:
            if isinstance(inlet, Stream) and inlet.sink is not self:
                raise ValueError("sink of given inlet must be this unit")
            else:
                outlet = self.outs[outlet]
            source.outs.replace(stream, inlet)

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
    
    def get_node(self):
        """Return unit node attributes for graphviz."""
        try: self._load_stream_links()
        except: pass
        if bst.preferences.minimal_nodes:
            return self._graphics.get_minimal_node(self)
        else:
            node = self._graphics.get_node_tailored_to_unit(self)
            node['tooltip'] = self._get_tooltip_string()
            return node
    
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
    
    @piping.ignore_docking_warnings
    def take_place_of(self, other, discard=False):
        """Replace inlets and outlets from this unit or system with that of 
        another unit or system."""
        self._ins[:] = other.ins
        self._outs[:] = other.outs
        if discard: bst.main_flowsheet.unit.discard(other)
    
    @piping.ignore_docking_warnings
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
        if discard: bst.main_flowsheet.unit.discard(self)
                    
    # Forward pipping
    def __sub__(self, other):
        """Source streams."""
        isa = isinstance
        if isa(other, (Unit, bst.System)):
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

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    def reset_cache(self, isdynamic=None):
        pass
    
    def _get_design_info(self): 
        return ()
        
    def _load_stream_links(self):
        options = self._stream_link_options
        if options:
            s_in = self._ins[0]
            s_out = self._outs[0]
            try:
                s_out.link_with(s_in, options.flow, options.phase, options.TP)
            except:
                pass
    
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
        self._specifications.append(bst.BoundedNumericalSpecification(f, *args, **kwargs))
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
        self._load_utility_cost()

    def _load_utility_cost(self):
        ins = self._ins._streams
        outs = self._outs._streams
        prices = bst.stream_utility_prices
        self._utility_cost = (
            sum([i.cost for i in self.heat_utilities]) 
            + self.power_utility.cost
            + sum([s.F_mass * prices[name] for name, index in self._inlet_utility_indices.items() if (s:=ins[index]).price == 0.])
            - sum([s.F_mass * prices[name] for name, index in self._outlet_utility_indices.items() if (s:=outs[index]).price == 0.])
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

    def results(self, with_units=True, include_utilities=True,
                include_total_cost=True, include_installed_cost=False,
                include_zeros=True, external_utilities=None, key_hook=None):
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
            if ' - ' in name:
                auxname, _ = name.split(' - ')
                auxname = auxname.replace(' ', '_').lower()
                for i in self.auxiliary_unit_names:
                    if auxname == i.lower():
                        N = int(parallel.get(i, N_default))
                        break
                else:
                    N = int(parallel.get(name, N_default))
            else:
                N = int(parallel.get(name, N_default))
            if N != 1: key = (*others, name + f' (x{N})')
            keys.append(key)
        parallel = self.parallel
        N_default = parallel.get('self', 1)
        keys = []; 
        vals = []; addval = vals.append
        stream_utility_prices = bst.stream_utility_prices
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
                for name, flow in self.get_inlet_utility_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (inlet)'
                        addkey((ID, 'Flow'))
                        addval(('kg/hr', flow))
                        addkey((ID, 'Cost'))
                        addval(('USD/hr', flow * stream_utility_prices[name]))
                for name, flow in self.get_outlet_utility_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (outlet)'
                        addkey((ID, 'Flow'))
                        addval(('kg/hr', flow))
                        addkey((ID, 'Cost'))
                        addval(('USD/hr', - flow * stream_utility_prices[name]))
                
            units = self._units
            Cost = self.purchase_costs
            for ki, vi in self.design_results.items():
                addkey(('Design', ki))
                addval((units.get(ki, ''), vi))
            for ki, vi, ui in self._get_design_info():
                addkey(('Design', ki))
                addval((ui, vi))
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
                for name, flow in self.get_inlet_utility_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (inlet)'
                        addkey((ID, 'Flow'))
                        addval(flow)
                        addkey((ID, 'Cost'))
                        addval(flow * stream_utility_prices[name])
                for name, flow in self.get_outlet_utility_flows().items():
                    if include_zeros or flow:
                        ID = name + ' (outlet)'
                        addkey((ID, 'Flow'))
                        addval(flow)
                        addkey((ID, 'Cost'))
                        addval(-flow * stream_utility_prices[name])
                            
            for ki, vi in self.design_results.items():
                addkey(('Design', ki))
                addval(vi)
            for ki, vi, ui in self._get_design_info():
                addkey(('Design', ki))
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

    @property
    def thermo(self) -> tmo.Thermo:
        """Thermodynamic property package."""
        return self._thermo
    @property
    def ins(self) -> Sequence[Stream]:
        """List of all inlet streams."""
        return self._ins    
    @property
    def outs(self) -> Sequence[Stream]:
        """List of all outlet streams."""
        return self._outs

    def get_available_chemicals(self):
        streams = [i for i in (self._ins + self._outs) if i]
        reaction_chemicals = sum([i.reaction_chemicals for i in self.__dict__.values() if hasattr(i, 'reaction_chemicals')], [])
        required_chemicals = set(sum([i.available_chemicals for i in streams], reaction_chemicals))
        return [i for i in self.chemicals if i in required_chemicals]

    def _add_upstream_neighbors_to_set(self, set, ends, facilities):
        """Add upsteam neighboring units to set."""
        for s in self._ins:
            u = s._source
            if u and (facilities or not isinstance(u, bst.Facility)) and not (ends and s in ends):
                set.add(u)

    def _add_downstream_neighbors_to_set(self, set, ends, facilities):
        """Add downstream neighboring units to set."""
        for s in self._outs:
            u = s._sink
            if u and (facilities or not isinstance(u, bst.Facility)) and not (ends and s in ends):
                set.add(u)

    def get_downstream_units(self, ends=None, facilities=True):
        """Return a set of all units downstream."""
        downstream_units = set()
        outer_periphery = set()
        self._add_downstream_neighbors_to_set(outer_periphery, ends, facilities)
        inner_periphery = None
        old_length = -1
        new_length = 0
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            downstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                unit._add_downstream_neighbors_to_set(outer_periphery, ends, facilities)
            new_length = len(downstream_units)
        return downstream_units
    
    def get_upstream_units(self, ends=None, facilities=True):
        """Return a set of all units upstream."""
        upstream_units = set()
        outer_periphery = set()
        self._add_upstream_neighbors_to_set(outer_periphery, ends, facilities)
        inner_periphery = None
        old_length = -1
        new_length = 0
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            upstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                unit._add_upstream_neighbors_to_set(outer_periphery, ends, facilities)
            new_length = len(upstream_units)
        return upstream_units
    
    def neighborhood(self, 
            radius: Optional[int]=1, 
            upstream: Optional[bool]=True,
            downstream: Optional[bool]=True, 
            ends: Optional[Stream]=None, 
            facilities: Optional[bool]=None
        ):
        """
        Return a set of all neighboring units within given radius.
        
        Parameters
        ----------
        radius : 
            Maximum number streams between neighbors.
        downstream : 
            Whether to include downstream operations.
        upstream : 
            Whether to include upstream operations.
        ends :
            Streams that mark the end of the neighborhood.
        facilities :
            Whether to include facilities.
        
        """
        radius -= 1
        neighborhood = set()
        if radius < 0: return neighborhood
        if upstream:self._add_upstream_neighbors_to_set(neighborhood, ends, facilities)
        if downstream: self._add_downstream_neighbors_to_set(neighborhood, ends, facilities)
        direct_neighborhood = neighborhood
        for i in range(radius):
            neighbors = set()
            for unit in direct_neighborhood:
                if upstream: unit._add_upstream_neighbors_to_set(neighbors, ends, facilities)
                if downstream: unit._add_downstream_neighbors_to_set(neighbors, ends, facilities)
            if neighbors == direct_neighborhood: break
            direct_neighborhood = neighbors
            neighborhood.update(direct_neighborhood)
        return neighborhood

    def diagram(self, radius: Optional[int]=0, upstream: Optional[bool]=True,
                downstream: Optional[bool]=True, file: Optional[str]=None, 
                format: Optional[str]=None, display: Optional[bool]=True,
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
        
        """
        if radius > 0:
            units = self.neighborhood(radius, upstream, downstream)
            units.add(self)
        else:
            units = [self]
        return bst.System(None, units).diagram(format=format, display=display, file=file, title='', **graph_attrs)
    
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
        return self._mol_in/self.F_mol_in
    @property
    def z_mol_out(self) -> NDArray[float]:
        """Molar fractions going in."""
        return self._mol_out/self.F_mol_out

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
        return sum(self.vol_in)

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
        return sum(self.vol_out)
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
    
    # Representation
    def _info(self, layout, T, P, flow, composition, N, IDs, sort, data):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        return info + repr_ins_and_outs(layout, self.ins, self.outs, T, P, flow, composition, N, IDs, sort, data)

    def show(self, layout=None, T=None, P=None, flow=None, composition=None, N=None, IDs=None, sort=None, data=True):
        """Prints information on unit."""
        print(self._info(layout, T, P, flow, composition, N, IDs, sort, data))
    
    def _ipython_display_(self):
        if bst.preferences.autodisplay: self.diagram()
        self.show()


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

del thermo_user, registered
