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
from thermosteam._graphics import UnitGraphics
from ._heat_utility import HeatUtility
from .utils import AbstractMethod, format_title, static, piping
from ._power_utility import PowerUtility
from .exceptions import UnitInheritanceError
from thermosteam.units_of_measure import convert
from copy import copy
import biosteam as bst
from thermosteam import Stream, AbstractUnit, ProcessSpecification
from numpy.typing import NDArray
from typing import Callable, Optional, TYPE_CHECKING, Sequence, Iterable
from thermosteam.base.sparse import SparseVector, sum_sparse_vectors
import numpy as np
import thermosteam as tmo
from thermosteam import EquationNode, VariableNode
if TYPE_CHECKING: 
    System = bst._recycle_system
    HXutility = bst.HXutility
    UtilityAgent = bst.UtilityAgent

streams = Optional[Sequence[Stream|str]]

__all__ = ('Unit',)

# %% Typing

# from typing import Collection, Union, Annotated
# streams = Union[Collection[Union[Stream, str, None]], Union[Stream, str, None]]
# stream = Union[Annotated[Union[Stream, str, None], 1], Union[Stream, str, None]]
# stream_sequence = Collection[Union[Stream, str, None]]

# %% Unit Operation

def phenomena_based_run(self):
    if not (self._recycle_system and self._system.algorithm == 'Phenomena based'):
        Unit.run(self)
        return
    ins = self.ins
    outs = self.outs
    Ts = [i.T for i in outs]
    Q = self.Hnet
    Unit.run(self)
    new = sum(
        [i.mol for i in outs],
        -sum([i.mol for i in ins], 0)
    )
    if hasattr(self, '_dmol'):
        old = self._dmol
        f = bst.PhasePartition.dmol_relaxation_factor
        self._dmol = dmol = f * old + (1 - f) * new
    else:
        self._dmol = dmol = new
    self._dmol[np.abs(dmol) < 1e-9] = 0.
    self._duty = self.Hnet
    Ts_new = [i.T for i in outs]
    if all([i == j for i, j in zip(Ts_new, Ts)]): # T constant
        self._energy_variable = None
    elif self.Hnet == Q: # Q constant    
        self._energy_variable = 'T'
    else:
        self._energy_variable = None

class Unit(AbstractUnit):
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
    Stream = Stream
    max_parallel_units = int(10e3)
    
    def __init_subclass__(cls,
                          isabstract=False,
                          new_graphics=True,
                          does_nothing=None,
                          default_phenomena=None):
        super().__init_subclass__()
        if does_nothing: return 
        dct = cls.__dict__
        if 'run' in dct:
            raise UnitInheritanceError(
                 "the 'run' method cannot be overridden; implement `_run` instead"
            )
        if default_phenomena is None:
            methods = (
                '_update_nonlinearities',
                '_update_energy_coefficient',
                '_update_variable',
                '_create_material_balance_equations',
                '_create_energy_balance_equations',
                '_create_bulk_balance_equations',
                '_energy_variable',
            )
            if all([getattr(cls, i) is getattr(Unit, i) for i in methods]):
                cls.run = phenomena_based_run
            else:
                for i in methods: 
                    if getattr(cls, i) is not getattr(Unit, i): continue
                    setattr(cls, i, AbstractMethod)
                cls.run = Unit.run
        elif default_phenomena:
            cls.run = phenomena_based_run
            for i in methods: setattr(cls, i, getattr(Unit, i))
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
        name = cls.__name__
        if hasattr(bst, 'units') and hasattr(bst, 'wastewater') and hasattr(bst, 'facilities'):
            # Add 3rd party unit to biosteam module for convenience
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

    #: **class-attribute** Whether to link inlet and outlet streams.
    _link_streams: bool = False

    #: **class-attribute** Lifetime of equipment. Defaults to lifetime of
    #: production venture. Use an integer to specify the lifetime for all
    #: items in the unit purchase costs. Use a dictionary to specify the 
    #: lifetime of each purchase cost item.
    _default_equipment_lifetime: int|dict[str, int] = {}

    #: [str] The energy variable for phenomena-based simulation.
    _energy_variable: str = None

    ### Abstract methods ###
    
    #: Create auxiliary components.
    _load_components = AbstractMethod
    
    #: Add design requirements to the :attr:`~Unit.design_results` dictionary.
    _design = AbstractMethod
    
    #: Add itemized purchase costs to the :attr:`~Unit.baseline_purchase_costs` dictionary.
    _cost = AbstractMethod

    #: Add embodied emissions (e.g., unit construction) in LCA
    _lca = AbstractMethod
    
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
        if self._energy_variable == 'T':
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
    
    Inlets = piping.Inlets
    Outlets = piping.Outlets
    
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
        
        super().__init__(ID, ins, outs, thermo, **kwargs)
    
        self._assert_compatible_property_package()
        
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
    
    def define_credit(self, name: str, stream: Stream):
        """
        Define an inlet or outlet stream as a fee/credit by name.
        
        Parameters
        ----------
        name : 
            Name of fee/credit, as defined in :meth:`settings.stream_prices <thermosteam._settings.ProcessSettings.stream_prices>`.
        stream :
            Inlet or outlet fee/credit stream.
        
        """
        if name not in bst.stream_prices:
            raise ValueError(f"price of '{name}' must be defined in settings.stream_prices")
        if stream._sink is self:
            self._inlet_cost_indices[name] = self._ins._streams.index(stream)
        elif stream._source is self:
            self._outlet_revenue_indices[name] = self._outs._streams.index(stream)
        else:
            raise ValueError(f"stream '{stream.ID}' must be connected to {repr(self)}")
    
    define_fee = define_credit
    
    def get_inlet_cost_flows(self):
        ins = self._ins._streams
        return {name: ins[index].F_mass for name, index in (self._inlet_utility_indices | self._inlet_cost_indices).items()}
    
    def get_outlet_revenue_flows(self):
        outs = self._outs._streams
        return {name: outs[index].F_mass for name, index in (self._outlet_utility_indices | self._outlet_revenue_indices).items()}
    
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
    
    def _get_design_info(self): 
        return ()
    
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
        self._inlet_cost = sum(
            [ins[index].F_mass * prices[name] for name, index in self._inlet_cost_indices.items()]
        )
        self._outlet_revenue = sum(
            [outs[index].F_mass * prices[name] for name, index in self._outlet_revenue_indices.items()]
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

    def _mass_and_energy_balance_specifications(self):
        return (self.line, ())

    def _mass_and_energy_balance_specifications_table(self):
        title, specs = self._mass_and_energy_balance_specifications()
        index = []
        values = []
        for name, value, units in specs:
            index.append(name)
            values.append([value, units])
        df = pd.DataFrame(
            values, index, ('Values', 'Units')
        )
        df.columns.name = f"{title} - {self.ID}"
        return df

    def results(self, with_units=True, include_utilities=True,
                include_total_cost=True, include_installed_cost=False,
                include_zeros=True, external_utilities=None, key_hook=None,
                basis=None):
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
                for ki, vi in self.design_results.items():
                    addkey(('Design', ki))
                    if ki in units:
                        ui = units[ki]
                        if ui != 'hr':
                            q = Quantity(vi, ui)
                            q = q.to_base_units()
                            vi = q.magnitude
                            ui = q.units
                    else:
                        ui = ''
                    addval((ui, vi))
                for ki, vi, ui in self._get_design_info():
                    addkey(('Design', ki))
                    q = Quantity(vi, ui)
                    q = q.to_base_units()
                    vi = q.magnitude
                    ui = q.units
                    addval((ui, vi))
            else:
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
