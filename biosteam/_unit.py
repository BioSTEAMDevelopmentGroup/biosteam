# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import pandas as pd
from warnings import warn
from graphviz import Digraph
from ._graphics import UnitGraphics, box_graphics
from thermosteam import Stream
from ._heat_utility import HeatUtility
from .utils import NotImplementedMethod, format_title, static
from .utils import piping
from ._power_utility import PowerUtility
from .digraph import finalize_digraph
from thermosteam.utils import thermo_user, registered
from thermosteam.units_of_measure import convert
from copy import copy
import biosteam as bst
import thermosteam as tmo

__all__ = ('Unit',)

_count = [0]
def count():
    _count[0] += 1
    print(_count)

# %% Inlet and outlet representation

def repr_ins_and_outs(layout, ins, outs, T, P, flow, composition, N, IDs, data):
    info = ''
    if ins:
        info += 'ins...\n'
        i = 0
        for stream in ins:
            unit = stream._source
            source_info = f'  from  {type(unit).__name__}-{unit}\n' if unit else '\n'
            if stream and data:
                stream_info = stream._info(layout, T, P, flow, composition, N, IDs)
                index = stream_info.index('\n')
                info += f'[{i}] {stream}' + source_info + stream_info[index+1:] + '\n'
            else:
                info += f'[{i}] {stream}' + source_info
            i += 1
    if outs:
        info += 'outs...\n'
        i = 0
        for stream in outs:
            unit = stream._sink
            sink_info = f'  to  {type(unit).__name__}-{unit}\n' if unit else '\n'
            if stream and data:
                stream_info = stream._info(layout, T, P, flow, composition, N, IDs)
                index = stream_info.index('\n')
                info += f'[{i}] {stream}' + sink_info + stream_info[index+1:] + '\n'
            else:
                info += f'[{i}] {stream}' + sink_info
            i += 1
    info = info.replace('\n ', '\n    ')
    return info[:-1]

# %% Path utilities

def find_path_segment(start_unit, end_unit):
    path_segment = fill_path_segment(start_unit, [], end_unit)
    if path_segment is None:
        raise ValueError(f"end unit {repr(end_unit)} not downstream from start unit {repr(start_unit)}")
    return path_segment

def fill_path_segment(start_unit, path, end_unit):
    if start_unit is end_unit:
        return path
    if start_unit not in path: 
        path.append(start_unit)
        first_outlet, *other_outlets = start_unit._outs
        for outlet in other_outlets:
            start_unit = outlet._sink
            if not start_unit: continue
            new_path = path.copy()
            path_segment = fill_path_segment(start_unit, new_path, end_unit)
            if path_segment: return path_segment
        start_unit = first_outlet._sink
        if not start_unit: return None
        path_segment = fill_path_segment(start_unit, path, end_unit)
        if path_segment: return path_segment
        

# %% Unit Operation

@thermo_user
@registered(ticket_name='U')
class Unit:
    """
    Abstract parent class for Unit objects. Child objects must contain
    `_run`, `_design` and `_cost` methods to estimate stream outputs of a
    Unit and find design and cost information.  

    **Abstract class methods**
    
    reset_cache()
        Reset unit operartion cache.
    _setup()
        Set stream conditions and constant data.
    _run()
        Run simulation and update output streams.
    _design()
        Add design requirements to the `design_results` dictionary.
    _cost()
        Add itemized purchse costs to the `baseline_purchase_costs` dictionary.

    **Abstract class attributes**
    
    **line='Unit'**
        [str] Name denoting the type of Unit class. Defaults to the class
        name of the first child class.
    **_F_BM_default** 
        dict[str, float] Default bare-module factors for each purchase cost item.
        Items in this dictionary are copied to the `F_BM` attribute during 
        initialization.
    **_units**
        [dict] Units of measure for `design_results` dictionary.
    **_N_ins=1**
        [int] Expected number of input streams.
    **_N_outs=2**
        [int] Expected number of output streams.
    **_ins_size_is_fixed=True**
        [bool] Whether the number of streams in ins is fixed.
    **_outs_size_is_fixed=True**
        [bool] Whether the number of streams in outs is fixed.
    **_N_heat_utilities=0**
        [int] Number of heat utilities created with each instance.
    **auxiliary_unit_names=()
        tuple[str] Name of attributes that are auxiliary units. These units
        will be accounted for in the purchase and installed equipment costs
        without having add these costs in the `purchase_costs` dictionary.
        Utility costs, however, are not automatically accounted for and must
        be hardcoded in the unit operation logic.
    **_default_equipment_lifetime=None**
        [int] or dict[str, int] Lifetime of equipment. Defaults to lifetime of
        production venture. Use an integer to specify the lifetime for all
        items in the unit purchase costs. Use a dictionary to specify the 
        lifetime of each purchase cost item.
    **_graphics**
        [biosteam.Graphics, abstract, optional] Settings for diagram
        representation. Defaults to a box with the same number of input
        and output edges as `_N_ins` and `_N_outs`.

    Parameters
    ----------
    ID='' : str, defaults to a unique ID
        A unique identification. If ID is None, unit will not be
        registered in flowsheet.
    ins=None : Iterable[:class:`~thermosteam.Stream`, or str], :class:`~thermosteam.Stream`, or str
        Inlet streams or IDs to initialize input streams.
        If empty, default IDs will be given. If None, defaults to missing streams.
    outs=() : Iterable[:class:`~thermosteam.Stream`, or str], :class:`~thermosteam.Stream`, or str
        Outlet streams or IDs to initialize output streams.
        If empty, default IDs will be given.
        If None, leave streams missing.
    thermo=None : :class:`~thermosteam.Thermo`
        Thermo object to initialize inlet and outlet streams. Defaults to
        `biosteam.settings.get_thermo()`.
    
    Attributes
    ----------
    ins : Inlets[:class:`~thermosteam.Stream`]
        Input streams.
    outs : Outlets[:class:`~thermosteam.Stream`]
        Output streams.
    power_utility : PowerUtility
        Electricity rate requirements are stored here (including auxiliary units).
    heat_utilities : tuple[:class:`~biosteam.HeatUtility`]
        Cooling and heating requirements are stored here (including auxiliary units).
    design_results : dict
        All design requirements (not including auxiliary units).
    purchase_costs : dict[str, float]
        Itemized purchase costs (including auxiliary units).
    thermo : Thermo
        The thermodynamic property package used by the unit.
    
    Examples
    --------
    :doc:`tutorial/Creating_a_Unit`
    
    :doc:`tutorial/Using_-pipe-_notation`
    
    :doc:`tutorial/Inheriting_from_Unit`
    
    :doc:`tutorial/Unit_decorators`
    
    """ 
    
    def __init_subclass__(cls,
                          isabstract=False,
                          new_graphics=True):
        dct = cls.__dict__
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
            elif 'mixer' in line: cls.ticket_name = 'M'
            elif 'pump' in line: cls.ticket_name = 'P'
            elif 'reactor' in line or 'digestion' in line or 'ferment' in line: cls.ticket_name = 'R'
            elif 'split' in line: cls.ticket_name = 'S'
            elif 'tank' in line: cls.ticket_name = 'T'
            elif 'junction' == line: cls.ticket_name = 'J'
            elif 'specification' in line: cls.ticket_name = 'PS'
            else: cls.ticket_name = 'U'
        if '_graphics' not in dct and new_graphics:
            # Set new graphics for specified line
            cls._graphics = UnitGraphics.box(cls._N_ins, cls._N_outs)
        if not isabstract:
            if cls is not Unit:
                if hasattr(cls, '_BM'): 
                    raise NotImplementedError(
                        'the `_BM` class attribute for bare-module factors is '
                        'deprecated; implement `_F_BM_default` instead'
                    )
                elif hasattr(cls, '_F_BM_defaults'):
                    raise NotImplementedError(
                        '`_F_BM_defaults` is incorrect; implement '
                        '`_F_BM_default` instead'
                    )
                elif not hasattr(cls, '_F_BM_default'):
                    cls._F_BM_default = {}
                
                if hasattr(cls, '_equipment_lifetime'):
                    raise NotImplementedError(
                        'the `_equipment_lifetime` class attribute is '
                        'deprecated; implement `_default_equipment_lifetime` instead'
                    )
                elif hasattr(cls, '_default_equipment_lifetimes'):
                    raise NotImplementedError(
                        '`_default_equipment_lifetimes` is incorrect; implement '
                        '`_default_equipment_lifetime` instead'
                    )
                elif not hasattr(cls, '_default_equipment_lifetime'): 
                    cls._default_equipment_lifetime = {}
            if not hasattr(cls, '_units'): cls._units = {}
            if not cls._run:
                if cls._N_ins == 1 and cls._N_outs == 1:
                    static(cls)
                else:
                    raise NotImplementedError(
                        "Unit subclass with multiple inlet or outlet streams "
                        "must implement a '_run' method unless the "
                        "'isabstract' keyword argument is True"
                    )
        if '__init__' in dct and '_stacklevel' not in dct:
            cls._stacklevel += 1
        
    ### Abstract Attributes ###
    
    # tuple[str] Name of attributes that are auxiliary units. These units
    # will be accounted for in the purchase and installed equipment costs
    # without having add these costs in the `purchase_costs` dictionary
    auxiliary_unit_names = ()
    
    # [int] Expected number of inlet streams
    _N_ins = 1  
    
    # [int] Expected number of outlet streams
    _N_outs = 1
    
    # [bool] Whether the number of streams in ins is fixed
    _ins_size_is_fixed = True
    
    # [bool] Whether the number of streams in outs is fixed
    _outs_size_is_fixed = True
    
    # [int] number of heat utilities
    _N_heat_utilities = 0
    
    # [StreamLinkOptions] Options for linking streams
    _stream_link_options = None
    
    # [biosteam Graphics] A Graphics object for diagram representation
    _graphics = box_graphics

    # [int] Used for piping warnings. Should be equal to 6 plus the number of
    # wrappers to Unit.__init__
    _stacklevel = 5
    
    # [str] The general type of unit, regardless of class
    line = 'Unit'

    ### Other defaults ###
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        self._register(ID)
        self._specification = None
        self._load_thermo(thermo)
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_utils()
        self._init_results()
        self._init_specification()
        self._assert_compatible_property_package()
    
    def _init_ins(self, ins):
        # Inlets[:class:`~thermosteam.Stream`] Input streams
        self._ins = piping.Inlets(self, self._N_ins, ins, self._thermo, 
                                  self._ins_size_is_fixed, self._stacklevel)
    
    def _init_outs(self, outs):
        # Outlets[:class:`~thermosteam.Stream`] Output streams
        self._outs = piping.Outlets(self, self._N_outs, outs, self._thermo,
                                    self._outs_size_is_fixed, self._stacklevel)
    
    def _init_utils(self):
        # tuple[HeatUtility] All heat utilities associated to unit
        self.heat_utilities = tuple([HeatUtility() for i in
                                     range(self._N_heat_utilities)])
        
        # [PowerUtility] Electric utility associated to unit
        self.power_utility = PowerUtility()
    
    def _init_results(self):
        try:
            #: [dict] All bare-module factors for each purchase cost.
            #: Defaults to values in the class attribute `_F_BM_default`.
            self.F_BM = self._F_BM_default.copy()
        except AttributeError:
            self.F_BM = {}
        
        #: [dict] All design factors for each purchase cost.
        self.F_D = {}
        
        #: [dict] All pressure factors for each purchase cost.
        self.F_P = {}
        
        #: [dict] All material factors for each purchase cost.
        self.F_M = {}
        
        # [dict] All design results.
        self.design_results = {}
        
        # [dict] All baseline purchase costs without accounting for design, 
        # pressure, and material factors.
        self.baseline_purchase_costs = {}
        
        # [dict] All purchase costs in USD.
        self.purchase_costs = {}
        
        # [dict] All installed costs accounting for bare module, design, 
        # pressure, and material factors.
        self.installed_costs = {}
        
        # dict[str: tuple(int, float)] Indices of additional utilities given by inlet streams.
        self.inlet_utility_indices = {}
        
        # dict[str: tuple(int, float)] Indices of additional utilities given by outlet streams.
        self.outlet_utility_indices = {}
        
        try:
            #: [int] or dict[str, int] Lifetime of equipment. Defaults to values in
            #: the class attribute `_default_equipment_lifetime`. Use an integer 
            #: to specify the lifetime for all items in the unit purchase costs.
            #: Use a dictionary to specify the lifetime of each purchase cost item.
            self.equipment_lifetime = copy(self._default_equipment_lifetime)
        except AttributeError:
            self.equipment_lifetime = {}
        
        #: [dict] Greenhouse gas emissions for use in BioSTEAM-LCA 
        #: (https://github.com/scyjth/biosteam_lca)
        self._GHGs = {}
    
    def _init_specification(self):
        #: list[Callable] All specification functions
        self._specification = []
        
        #: [bool] Whether to run mass and energy balance after calling
        #: specification functions
        self.run_after_specification = False 
    
    def _reset_thermo(self, thermo):
        for i in (self._ins._streams + self._outs._streams):
            try:
                if i: i._reset_thermo(thermo)
            except:
                raise RuntimeError(f'failed to reset {repr(self)}.thermo')
        if thermo is self.thermo: return
        self._load_thermo(thermo)
        chemicals = thermo.chemicals
        dcts = [self.__dict__]
        if hasattr(self, 'components') and isinstance(self.components, dict):
            dcts.append(self.components)
        for dct in dcts:
            for i, j in dct.items():
                if isinstance(j, tmo.ReactionSystem):
                    for rxn in j._reactions:
                        if hasattr(rxn, 'reset_chemicals') and rxn.chemicals is not chemicals:
                            rxn.reset_chemicals(chemicals)
                elif hasattr(j, 'reset_chemicals') and j.chemicals is not chemicals:
                    j.reset_chemicals(chemicals)
                elif hasattr(j, '_reset_thermo') and j.thermo is not thermo:
                    j._reset_thermo(thermo)
    
    
    def define_utility(self, name, stream):
        if name not in bst.stream_utility_prices:
            raise ValueError(f"price of '{name}' must be defined in biosteam.stream_utility_prices")
        if stream._sink is self:
            self.inlet_utility_indices[name] = self._ins._streams.index(stream)
        elif stream._source is self:
            self.outlet_utility_indices[name] = self._outs._streams.index(stream)
        else:
            raise ValueError(f"stream '{stream.ID}' must be connected to {repr(self)}")
            
    def get_inlet_utility_flows(self):
        ins = self._ins._streams
        return {name: ins[index].F_mass for name, index in self.inlet_utility_indices.items()}
    
    def get_outlet_utility_flows(self):
        outs = self._outs._streams
        return {name: outs[index].F_mass for name, index in self.outlet_utility_indices.items()}
    
    def get_design_and_capital(self):
        return UnitDesignAndCapital(
            self, self.F_BM.copy(), self.F_D.copy(), self.F_P.copy(), self.F_M.copy(), 
            self.design_results.copy(), self.baseline_purchase_costs.copy(),
            self.purchase_costs.copy(), self.installed_costs.copy(),
        )
    
    def get_agile_design_and_capital(self, design_and_capital: list):
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
                self._load_capital_costs()
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
        for name, Cp in baseline_purchase_costs.items():
            F = F_D.get(name, 1.) * F_P.get(name, 1.) * F_M.get(name, 1.)
            installed_cost = Cp * (F_BM.get(name, 1.) + F - 1.)
            purchase_cost = Cp * F
            if installed_cost > installed_costs[name]:
                installed_costs[name] = installed_cost
            if purchase_cost > purchase_costs[name]:
                purchase_costs[name] = purchase_cost
    
    def _assert_compatible_property_package(self):
        chemicals = self.chemicals
        streams = self._ins + self._outs
        assert all([s.chemicals is chemicals for s in streams if s]), (
            "unit operation chemicals are incompatible with inlet and outlet streams; "
            "try using the `thermo` keyword argument to initialize the unit operation "
            "with a compatible thermodynamic property package"
        )
    
    def _load_auxiliary_capital_costs(self):
        for i in self.auxiliary_units: i._load_capital_costs()
        baseline_purchase_costs = self.baseline_purchase_costs
        purchase_costs = self.purchase_costs
        installed_costs = self.installed_costs
        F_BM = self.F_BM
        F_D = self.F_D
        F_P = self.F_P
        F_M = self.F_M
        for name in self.auxiliary_unit_names:
            unit = getattr(self, name)
            if not unit: continue
            F_BM_auxiliary = unit.F_BM
            F_D_auxiliary = unit.F_D
            F_P_auxiliary = unit.F_P
            F_M_auxiliary = unit.F_M
            bpc_auxiliary = unit.baseline_purchase_costs
            pc_auxiliary = unit.purchase_costs
            ic_auxiliary = unit.installed_costs
            for i in unit.baseline_purchase_costs:
                j = ' - '.join([name.capitalize(), i])
                if j in baseline_purchase_costs: 
                    raise RuntimeError(
                        f"'{j}' already in `baseline_purchase_cost` "
                        f"dictionary of {repr(self)}; try using a different key"
                    )
                else:
                    F_D[j] = fd = F_D_auxiliary.get(i, 1.)
                    F_P[j] = fp = F_P_auxiliary.get(i, 1.)
                    F_M[j] = fm = F_M_auxiliary.get(i, 1.)
                    baseline_purchase_costs[j] = Cpb = bpc_auxiliary[i]
                    purchase_costs[j] = pc_auxiliary[i]
                    installed_costs[j] = Cbm = ic_auxiliary[i]
                    try:
                        F_BM[j] = F_BM_auxiliary[i]
                    except KeyError:
                        # Assume costs already added elsewhere using another method.
                        # Calculate BM as an estimate.
                        F_BM[j] = Cbm / Cpb + 1 - fd * fp * fm
    
    def _load_capital_costs(self):
        r"""
        Calculate and save free on board (f.o.b.) purchase costs and
        installed equipment costs (i.e. bare-module cost) for each item in the 
        `baseline_purchase_costs` dictionary and in auxiliary units.
        
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
        each equipment should be stored in the `F_BM`, `F_D`, `F_P`, and 
        `F_M` dictionaries.
        
        Warning
        -------
        If an item is listed in the `purchase_costs` dictionary but not in the
        `baseline_purchase_costs` dictionary, the baseline purchase cost is 
        assumed to be the same as the purchase cost.
        
        References
        ----------
        .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
        
        """
        F_BM = self.F_BM
        F_D = self.F_D
        F_P = self.F_P
        F_M = self.F_M
        self._load_auxiliary_capital_costs()
        baseline_purchase_costs = self.baseline_purchase_costs
        purchase_costs = self.purchase_costs
        installed_costs = self.installed_costs
        for i in purchase_costs:
            if i not in baseline_purchase_costs:
                warning = RuntimeWarning(
                    "adding items to the `purchase_costs` dictionary is "
                    "deprecated; add items to `baseline_purchase_costs` "
                    "dictionary instead"
                 )
                warn(warning)
                baseline_purchase_costs[i] = purchase_costs[i]
        for name, Cpb in baseline_purchase_costs.items(): 
            if name in installed_costs and name in purchase_costs: 
                continue # Assume costs already added elsewhere using another method
            F = F_D.get(name, 1.) * F_P.get(name, 1.) * F_M.get(name, 1.)
            try:
                installed_costs[name] = Cpb * (F_BM[name] + F - 1.)
            except KeyError:
                warning = RuntimeWarning(
                   f"the purchase cost item, '{name}', has "
                    "no defined bare-module factor in the "
                  f"'{type(self).__name__}.F_BM' dictionary; "
                   "bare-module factor now has a default value of 1"
                 )
                warn(warning)
                F_BM[name] = 1.
                installed_costs[name] = purchase_costs[name] = Cpb * F
            else:
                purchase_costs[name] = Cpb * F
    
    def _setup(self):
        self.baseline_purchase_costs.clear()
        self.purchase_costs.clear()
        self.installed_costs.clear()
    
    def materialize_connections(self):
        for s in self._ins + self._outs: 
            if not s: s.materialize_connection()
    
    @property
    def owner(self):
        owner = getattr(self, '_owner', None)
        if owner is None:
            return self
        else:
            return owner.owner
    @owner.setter
    def owner(self, owner):
        if owner is self: return
        self._owner = owner
    
    def disconnect(self):
        self._ins[:] = ()
        self._outs[:] = ()
    
    def get_node(self):
        """Return unit node attributes for graphviz."""
        try: self._load_stream_links()
        except: pass
        if bst.MINIMAL_UNIT_DIAGRAMS:
            return self._graphics.get_minimal_node(self)
        else:
            return self._graphics.get_node_tailored_to_unit(self)
    
    def get_design_result(self, key, units):
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
        int_types = (int, np.int)
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
        int_types = (int, np.int)
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
    
    # Abstract methods
    reset_cache = NotImplementedMethod
    _load_components = NotImplementedMethod
    _run = NotImplementedMethod
    _design = NotImplementedMethod
    _cost = NotImplementedMethod
    
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
    
    def add_specification(self, specification=None, run=None):
        if not specification: return lambda specification: self.add_specification(specification, run)
        if not callable(specification): raise ValueError('specification must be callable')
        self.specification.append(specification)
        if run is not None: self.run_after_specification = run
        return specification
    
    def run(self):
        """Run mass and energy balance."""
        specification = self._specification
        if specification:
            for i in specification: i()
            if self.run_after_specification: self._run()
        else:
            self._run()
            
    def path_until(self, unit):
        """Return a list of units starting from this one until the end unit (not inclusive) by
        moving along unit connections."""
        return find_path_segment(self, unit)
            
    def _reevaluate(self):
        """Reevaluate design and costs."""
        self.baseline_purchase_costs.clear()
        self.purchase_costs.clear()
        self.installed_costs.clear()
        self._summary()
    
    def _summary(self):
        """Calculate all results from unit run."""
        if not (self._design or self._cost): return
        self._design()
        self._cost()
        self._load_capital_costs()
        ins = self._ins._streams
        outs = self._outs._streams
        prices = bst.stream_utility_prices
        self._utility_cost = (
            sum([i.cost for i in self.heat_utilities]) 
            + self.power_utility.cost
            + sum([ins[index].F_mass * prices[name] for name, index in self.inlet_utility_indices.items()])
            - sum([outs[index].F_mass * prices[name] for name, index in self.outlet_utility_indices.items()])
        )
    
    @property
    def specification(self):
        """Process specification."""
        return self._specification
    @specification.setter
    def specification(self, specification):
        if specification:
            if callable(specification):
                self._specification = specification = [specification]
            elif (isinstance(specification, list)
                  and all([callable(i) for i in specification])):
                self._specification = specification
            else:
                raise AttributeError(
                    "specification must be callable or a list of callables; "
                   f"not a '{type(specification).__name__}'"
                )
    
    @property
    def baseline_purchase_cost(self):
        """Total baseline purchase cost, without accounting for design ,
        pressure, and material factors [USD]."""
        return sum(self.baseline_purchase_costs.values())
    
    @property
    def purchase_cost(self):
        """Total purchase cost [USD]."""
        return sum(self.purchase_costs.values())
    
    @property
    def installed_cost(self):
        """Total installed equipment cost [USD]."""
        return sum(self.installed_costs.values())
    
    @property
    def utility_cost(self):
        """Total utility cost [USD/hr]."""
        try:
            return self._utility_cost
        except:
            ins = self._ins._streams
            outs = self._outs._streams
            prices = bst.stream_utility_prices
            self._utility_cost = (
                sum([i.cost for i in self.heat_utilities]) 
                + self.power_utility.cost
                + sum([ins[index].F_mass * prices[name] for name, index in self.inlet_utility_indices.items()])
                - sum([outs[index].F_mass * prices[name] for name, index in self.outlet_utility_indices.items()])
            )
            return self._utility_cost

    @property
    def auxiliary_units(self):
        """tuple[Unit] All associated auxiliary units."""
        getfield = getattr
        return tuple([getfield(self, i) for i in self.auxiliary_unit_names])

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

    def simulate(self):
        """
        Run rigorous simulation and determine all design requirements.
        
        """
        self._setup()
        self._load_stream_links()
        self.run()
        self._summary()

    def results(self, with_units=True, include_utilities=True,
                include_total_cost=True, include_installed_cost=False,
                include_zeros=True, external_utilities=(), key_hook=None):
        """
        Return key results from simulation as a DataFrame if `with_units`
        is True or as a Series otherwise.
        """
        # TODO: Divide this into functions
        def addkey(key):
            if key_hook: key = key_hook(key)
            keys.append(key)
        
        keys = []; 
        vals = []; addval = vals.append
        stream_utility_prices = bst.stream_utility_prices
        if with_units:
            if include_utilities:
                power_utility = self.power_utility
                if power_utility:
                    addkey(('Power', 'Rate'))
                    addval(('kW', power_utility.rate))
                    if include_zeros or power_utility.cost:
                        addkey(('Power', 'Cost'))
                        addval(('USD/hr', power_utility.cost))
                for heat_utility in HeatUtility.sum_by_agent(self.heat_utilities + external_utilities):
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
                addkey(('Purchase cost', ki))
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
            if self._GHGs:
                a, b = self._totalGHG
                GHG_units =  self._GHG_units
                for ko, vo in self._GHGs.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval((GHG_units.get(ko, ''), vi))
                a_key, b_key = GHG_units.keys()
                a_unit, b_unit = GHG_units.values()
                addkey(('Total ' + a_key, ''))
                addval((a_unit, a))
                addkey(('Total ' + b_key, ''))
                addval((b_unit, b))
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
                    addkey(('Power', 'Rate'))
                    addval(power_utility.rate)
                    if include_zeros or power_utility.cost:
                        addkey(('Power', 'Cost'))
                        addval(power_utility.cost)
                for heat_utility in HeatUtility.sum_by_agent(self.heat_utilities + external_utilities):
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
                addkey(('Purchase cost', ki))
                addval(vi)
            if self._GHGs:
                GHG_units = self._GHG_units
                for ko, vo in self._GHGs.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval(vi)
                a, b = self._totalGHG
                a_key, b_key = GHG_units.keys()
                addkey(('Total ' + a_key, ''))
                addval(a)
                addkey(('Total ' + b_key, ''))
                addval(b)
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
    def thermo(self):
        """Thermodynamic property package"""
        return self._thermo
    @property
    def ins(self):
        """All input streams."""
        return self._ins    
    @property
    def outs(self):
        """All output streams."""
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
    
    def neighborhood(self, radius=1, upstream=True, downstream=True, ends=None, facilities=None):
        """
        Return a set of all neighboring units within given radius.
        
        Parameters
        ----------
        radius : int
                 Maxium number streams between neighbors.
        downstream=True : bool, optional
            Whether to include downstream operations
        upstream=True : bool, optional
            Whether to include upstream operations
        
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

    def get_digraph(self, format='png', **graph_attrs):
        ins = self.ins
        outs = self.outs
        graphics = self._graphics

        # Make a Digraph handle
        f = Digraph(name='unit', filename='unit', format=format)
        f.attr('graph', ratio='0.5', outputorder='edgesfirst',
               nodesep='1.1', ranksep='0.8', maxiter='1000')  # Specifications
        f.attr(rankdir='LR', **graph_attrs)  # Left to right

        # If many streams, keep streams close
        if (len(ins) >= 3) or (len(outs) >= 3):
            f.attr('graph', nodesep='0.4')

        # Initialize node arguments based on unit and make node
        node = graphics.get_node_tailored_to_unit(self)
        f.node(**node)

        # Set stream node attributes
        f.attr('node', shape='rarrow', fillcolor='#79dae8',
               style='filled', orientation='0', width='0.6',
               height='0.6', color='black', peripheries='1')

        # Make nodes and edges for input streams
        di = 0  # Destination position of stream
        for stream in ins:
            if not stream: continue
            f.node(stream.ID)
            edge_in = graphics.edge_in
            if di >= len(edge_in): di = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(stream.ID, node['name'])
            di += 1

        # Make nodes and edges for output streams
        oi = 0  # Origin position of stream
        for stream in outs:
            if not stream: continue
            f.node(stream.ID) 
            edge_out = graphics.edge_out  
            if oi >= len(edge_out): oi = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(node['name'], stream.ID)
            oi += 1
        return f

    def diagram(self, radius=0, upstream=True, downstream=True, 
                file=None, format='png', display=True, **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram
        of the unit and all neighboring units within given radius.
        
        Parameters
        ----------
        radius : int
                 Maxium number streams between neighbors.
        downstream=True : bool, optional
            Whether to show downstream operations
        upstream=True : bool, optional
            Whether to show upstream operations
        file : Must be one of the following:
            * [str] File name to save diagram.
            * [None] Display diagram in console.
        format : str
                 Format of file.
        display : bool, optional
            Whether to display diagram in console or to return the graphviz 
            object.
        
        """
        if radius > 0:
            neighborhood = self.neighborhood(radius, upstream, downstream)
            neighborhood.add(self)
            sys = bst.System('', neighborhood)
            return sys.diagram('thorough', file, format, **graph_attrs)
        f = self.get_digraph(format, **graph_attrs)
        if display or file: 
            finalize_digraph(f, file, format)
        else:
            return f
    
    ### Net input and output flows ###
    
    # Molar flow rates
    @property
    def mol_in(self):
        """Molar flows going in [kmol/hr]."""
        return sum([s.mol for s in self._ins if s])
    @property
    def mol_out(self):
        """Molar flows going out [kmol/hr]."""
        return sum([s.mol for s in self._outs if s])

    @property
    def z_mol_in(self):
        """Molar fractions going in [kmol/hr]."""
        return self._mol_in/self.F_mol_in
    @property
    def z_mol_out(self):
        """Molar fractions going in."""
        return self._mol_out/self.F_mol_out

    @property
    def F_mol_in(self):
        """Net molar flow going in [kmol/hr]."""
        return sum([s.F_mol for s in self._ins if s])
    @property
    def F_mol_out(self):
        """Net molar flow going out [kmol/hr]."""
        return sum([s.F_mol for s in self._outs if s])

    # Mass flow rates
    @property
    def mass_in(self):
        """Mass flows going in [kg/hr]."""
        return sum([s.mol for s in self._ins if s]) * self._thermo.chemicals.MW
    @property
    def mass_out(self):
        """Mass flows going out [kg/hr]."""
        return sum([s.mol for s in self._outs if s]) * self._thermo.chemicals.MW

    @property
    def z_mass_in(self):
        """Mass fractions going in."""
        return self.mass_in/self.F_mass_in
    @property
    def z_mass_out(self):
        """Mass fractions going out."""
        return self.mass_out/self.F_mass_out

    @property
    def F_mass_in(self):
        """Net mass flow going in [kg/hr]."""
        return self.mass_in.sum()
    @property
    def F_mass_out(self):
        """Net mass flow going out [kg/hr]."""
        return self.mass_out.sum()

    # Volumetric flow rates
    @property
    def vol_in(self):
        """Volumetric flows going in [m3/hr]."""
        return sum([s.vol for s in self._ins if s])
    @property
    def F_vol_in(self):
        """Net volumetric flow going in [m3/hr]."""
        return sum(self.vol_in)

    @property
    def z_vol_in(self):
        """Volumetric fractions going in."""
        return self.vol_in/self.F_vol_in
    @property
    def vol_out(self):
        """Volumetric flows going out [m3/hr]."""
        return sum([s.vol for s in self._outs if s])

    @property
    def F_vol_out(self):
        """Net volumetric flow going out [m3/hr]."""
        return sum(self.vol_out)
    @property
    def z_vol_out(self):
        """Volumetric fractions going out."""
        return self.vol_out/self.F_vol_out

    # Enthalpy flow rates
    @property
    def H_in(self):
        """Enthalpy flow going in [kJ/hr]."""
        return sum([s.H for s in self._ins if s])

    @property
    def H_out(self):
        """Enthalpy flow going out [kJ/hr]."""
        return sum([s.H for s in self._outs if s])

    @property
    def Hf_in(self):
        """Enthalpy of formation flow going in [kJ/hr]."""
        return sum([s.Hf for s in self._ins if s])

    @property
    def Hf_out(self):
        """Enthalpy of formation flow going out [kJ/hr]."""
        return sum([s.Hf for s in self._outs if s])

    @property
    def Hnet(self):
        """Net enthalpy flow, including enthalpies of formation [kJ/hr]."""
        return self.H_out - self.H_in + self.Hf_out - self.Hf_in
    
    # Representation
    def _info(self, layout, T, P, flow, composition, N, IDs, data):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        return info + repr_ins_and_outs(layout, self.ins, self.outs, T, P, flow, composition, N, IDs, data)

    def show(self, layout=None, T=None, P=None, flow=None, composition=None, N=None, IDs=None, data=True):
        """Prints information on unit."""
        print(self._info(layout, T, P, flow, composition, N, IDs, data))
    
    def _ipython_display_(self):
        if bst.ALWAYS_DISPLAY_DIAGRAMS: self.diagram()
        self.show()

    def __repr__(self):
        if self.ID:
            return f'<{type(self).__name__}: {self.ID}>'
        else:
            return f'<{type(self).__name__}>'


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