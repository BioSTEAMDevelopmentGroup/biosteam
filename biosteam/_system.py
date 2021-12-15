# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                          Sarang Bhagwat <sarangb2@illinois.edu>,
#                          Joy Zhang <joycheung1994@gmail.com>,
#                          Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
from .digraph import (digraph_from_units_and_streams,
                      digraph_from_system,
                      minimal_digraph,
                      surface_digraph,
                      finalize_digraph)
from thermosteam import functional as fn
from thermosteam import Stream, MultiStream
from thermosteam.utils import registered
from .exceptions import try_method_with_object_stamp
from ._network import Network
from ._facility import Facility
from ._unit import Unit, repr_ins_and_outs
from .utils import repr_items, ignore_docking_warnings
from .report import save_report
from .exceptions import InfeasibleRegion
from .utils import StreamPorts, OutletPort, colors
from .process_tools import utils
# from .utils import NotImplementedMethod
from collections.abc import Iterable
from warnings import warn
from inspect import signature
from thermosteam.utils import repr_kwargs
import biosteam as bst
import numpy as np
from scipy.integrate import solve_ivp, odeint

__all__ = ('System', 'AgileSystem', 'MockSystem',
           'AgileSystem', 'OperationModeResults',
           'mark_disjunction', 'unmark_disjunction')

def _reformat(name):
    name = name.replace('_', ' ')
    if name.islower(): name= name.capitalize()
    return name

# %% LCA

class ProcessImpactItem:
    __slots__ = ('name', 'basis', 'inventory', 'CF')
    
    def __init__(self, name, basis, inventory, CF):
        self.name = name
        self.basis = basis
        self.inventory = inventory
        self.CF = CF
        
    def impact(self):
        return self.inventory() * self.CF


# %% Customization to system creation

disjunctions = []

def mark_disjunction(stream):
    port = OutletPort.from_outlet(stream)
    if port not in disjunctions:
        disjunctions.append(port)

def unmark_disjunction(stream):
    port = OutletPort.from_outlet(stream)
    if port in disjunctions:
        disjunctions.remove(port)


# %% Functions for creating deterministic systems

def facilities_from_units(units):
    isa = isinstance
    return [i for i in units if isa(i, Facility)]

def find_blowdown_recycle(facilities):
    isa = isinstance
    for i in facilities:
        if isa(i, bst.BlowdownMixer): return i.outs[0]

# %% Functions for taking care of numerical specifications within a system path

def converge_system_in_path(system):
    specification = system._specification
    if specification:
        method = specification
    else:
        method = system._converge
    try_method_with_object_stamp(system, method)

def simulate_unit_in_path(unit):
    try_method_with_object_stamp(unit, unit.simulate)


# %% Debugging and exception handling

def raise_recycle_type_error(recycle):
    raise ValueError(
       f"invalid recycle of type '{type(recycle).__name__}' encountered; "
        "recycle must be either a Stream object, a tuple of Stream objects, or None"
    )

def print_exception_in_debugger(self, func, e):
    print(f"{colors.exception(type(e).__name__+ ':')} {e}")
    try: self.show()
    except: pass

def update_locals_with_flowsheet(lcs):
    lcs.update(bst.main_flowsheet.to_dict())
    lcs.update(bst.__dict__)

def _method_debug(self, f):
    """Method decorator for debugging system."""
    def g(*args, **kwargs):
        try:
            f(*args, **kwargs)
        except Exception as e:
            print_exception_in_debugger(self, f, e)
            update_locals_with_flowsheet(locals())
            # All systems, units, streams, and flowsheets are available as
            # local variables. Although this debugging method is meant
            # for internal development, please feel free to give it a shot.
            breakpoint()
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    g._original = f
    return g

def _method_profile(self, f):
    self._total_excecution_time_ = 0.
    t = bst.utils.TicToc()
    def g():
        t.tic()
        f()
        self._total_excecution_time_ += t.elapsed_time
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    g._original = f
    return g

# %% Converging recycle systems

class MockSystem:
    """
    Create a MockSystem object with inlets and outlets just like System
    objects, but without implementing any of the convergence methods nor
    path related attributes.

    Parameters
    ----------
    units : Iterable[:class:`~biosteam.Unit`], optional
        Unit operations in mock system.

    Notes
    -----
    This object is used to prevent the creation of unneeded systems for less
    computational effort.

    """
    __slots__ = ('units',
                 'flowsheet',
                 '_ins',
                 '_outs',
                 '_irrelevant_units')

    def __init__(self, units=()):
        self.units = units or list(units)
        self._load_flowsheet()

    def _load_flowsheet(self):
        self.flowsheet = flowsheet_module.main_flowsheet.get_flowsheet()

    @property
    def ins(self):
        """StreamPorts[:class:`~InletPort`] All inlets to the system."""
        try:
            return self._ins
        except:
            inlets = bst.utils.feeds_from_units(self.units)
            self._ins = ins = StreamPorts.from_inlets(inlets, sort=True)
            return ins
    @property
    def outs(self):
        """StreamPorts[:class:`~OutletPort`] All outlets to the system."""
        try:
            return self._outs
        except:
            outlets = bst.utils.products_from_units(self.units)
            self._outs = outs = StreamPorts.from_outlets(outlets, sort=True)
            return outs

    def load_inlet_ports(self, inlets, optional=()):
        """Load inlet ports to system."""
        all_inlets = bst.utils.feeds_from_units(self.units)
        inlets = list(inlets)
        for i in inlets:
            if i not in all_inlets:
                if i in optional:
                    inlets.remove(i)
                else:
                    raise ValueError(f'{i} is not an inlet')
        self._ins = StreamPorts.from_inlets(inlets)

    def load_outlet_ports(self, outlets, optional=()):
        """Load outlet ports to system."""
        all_outlets = bst.utils.products_from_units(self.units)
        outlets = list(outlets)
        for i in outlets:
            if i not in all_outlets:
                if i in optional:
                    outlets.remove(i)
                else:
                    raise ValueError(f'{i} is not an outlet')
        self._outs = StreamPorts.from_outlets(outlets)

    def __enter__(self):
        if self.units:
            raise RuntimeError("only empty mock systems can enter `with` statement")
        unit_registry = self.flowsheet.unit
        self._irrelevant_units = set(unit_registry)
        unit_registry._open_dump(self)
        return self

    def __exit__(self, type, exception, traceback):
        irrelevant_units = self._irrelevant_units
        del self._irrelevant_units
        if self.units:
            raise RuntimeError('mock system was modified before exiting `with` statement')
        unit_registry = self.flowsheet.unit
        dump = unit_registry._close_dump(self)
        self.units = [i for i in dump if i not in irrelevant_units]
        if exception: raise exception

    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__
    __pow__ = __sub__
    __rpow__ = __rsub__

    def show(self):
        ins = repr_items('    ins=', self.ins._ports, brackets='[]')
        outs = repr_items('    outs=', self.outs._ports, brackets='[]')
        units = repr_items('    units=', self.units, brackets='[]')
        args = ',\n'.join([ins, outs, units])
        print(f"{type(self).__name__}(\n{args}\n)")

    _ipython_display_ = show

    def __repr__(self):
        return f"{type(self).__name__}(ins={self.ins}, outs={self.outs})"


@registered('SYS')
class System:
    """
    Create a System object that can iteratively run each element in a path
    of BioSTREAM objects until the recycle stream is converged. A path can
    have function, Unit and/or System objects. When the path contains an
    inner System object, it converges/solves it in each loop/iteration.

    Parameters
    ----------
    ID : str
         A unique identification. If ID is None, instance will not be
         registered in flowsheet.
    path : tuple[:class:`~biosteam.Unit`, function and/or :class:`~biosteam.System`], optional
        A path that is run element by element until the recycle converges.
    recycle=None : :class:`~thermosteam.Stream` or tuple[:class:`~thermosteam.Stream`], optional
        A tear stream for the recycle loop.
    facilities=() : tuple[:class:`~biosteam.Unit`, function, and/or :class:`~biosteam.System`], optional
        Offsite facilities that are simulated only after
        completing the path simulation.
    facility_recycle : :class:`~thermosteam.Stream`, optional
        Recycle stream between facilities and system path.
    N_runs : int, optional
        Number of iterations to run system. This parameter is applicable
        only to systems with no recycle loop.
    operating_hours : float, optional
        Number of operating hours in a year. This parameter is used to
        compute convinience properties such as utility cost and material cost
        on a per year basis.
    lang_factor : float, optional
        Lang factor for getting fixed capital investment from
        total purchase cost. If no lang factor, installed equipment costs are
        estimated using bare module factors.

    """
    __slots__ = (
        '_ID',
        '_path',
        '_facilities',
        '_facility_loop',
        '_recycle',
        '_N_runs',
        '_specification',
        '_mol_error',
        '_T_error',
        '_rmol_error',
        '_rT_error',
        '_iter',
        '_ins',
        '_outs',
        'maxiter',
        'molar_tolerance',
        'relative_molar_tolerance',
        'temperature_tolerance',
        'relative_temperature_tolerance',
        'operating_hours',
        'flowsheet',
        'lang_factor',
        'process_impact_items',
        '_stabilized',
        '_connections',
        '_irrelevant_units',
        '_converge_method',
        '_TEA',
        '_LCA',
        '_subsystems',
        '_units',
        '_unit_path',
        '_cost_units',
        '_streams',
        '_feeds',
        '_products',
        '_state',
    )

    take_place_of = Unit.take_place_of
    replace_with = Unit.replace_with

    ### Class attributes ###

    #: [int] Default maximum number of iterations
    default_maxiter = 200

    #: [float] Default molar tolerance for each component (kmol/hr)
    default_molar_tolerance = 1.

    #: [float] Default relative molar tolerance for each component
    default_relative_molar_tolerance = 0.01

    #: [float] Default temperature tolerance (K)
    default_temperature_tolerance = 0.10

    #: [float] Default relative temperature tolerance
    default_relative_temperature_tolerance = 0.001

    #: [str] Default convergence method.
    default_converge_method = 'Aitken'

    # [bool] Whether to use stabilized convergence algorithm.
    default_stabilized_convergence = False

    #: [bool] Whether to raise a RuntimeError when system doesn't converge
    strict_convergence = True

    @classmethod
    def from_feedstock(cls, ID, feedstock, feeds=None, facilities=(),
                       ends=None, facility_recycle=None, operating_hours=None,
                       lang_factor=None):
        """
        Create a System object from a feedstock.

        Parameters
        ----------
        ID : str
            Name of system.
        feedstock : :class:`~thermosteam.Stream`
            Main feedstock of the process.
        feeds : Iterable[:class:`~thermosteam.Stream`]
            Additional feeds to the process.
        facilities : Iterable[Facility]
            Offsite facilities that are simulated only after
            completing the path simulation.
        ends : Iterable[:class:`~thermosteam.Stream`]
            Streams that not products, but are ultimately specified through
            process requirements and not by its unit source.
        facility_recycle : [:class:`~thermosteam.Stream`], optional
            Recycle stream between facilities and system path.
        operating_hours : float, optional
            Number of operating hours in a year. This parameter is used to
            compute convinience properties such as utility cost and material cost
            on a per year basis.
        lang_factor : float, optional
            Lang factor for getting fixed capital investment from
            total purchase cost. If no lang factor, installed equipment costs are
            estimated using bare module factors.

        """
        network = Network.from_feedstock(feedstock, feeds, ends)
        return cls.from_network(ID, network, facilities,
                                facility_recycle, operating_hours,
                                lang_factor)

    @classmethod
    def from_units(cls, ID="", units=None, feeds=None, ends=None,
                   facility_recycle=None, operating_hours=None,
                   lang_factor=None):
        """
        Create a System object from all units and streams defined in the flowsheet.

        Parameters
        ----------
        ID : str, optional
            Name of system.
        units : Iterable[:class:`biosteam.Unit`], optional
            Unit operations to be included.
        feeds : Iterable[:class:`~thermosteam.Stream`], optional
            All feeds to the system. Specify this argument if only a section
            of the complete system is wanted as it may disregard some units.
        ends : Iterable[:class:`~thermosteam.Stream`], optional
            End streams of the system which are not products. Specify this
            argument if only a section of the complete system is wanted, or if
            recycle streams should be ignored.
        facility_recycle : :class:`~thermosteam.Stream`, optional
            Recycle stream between facilities and system path. This argument
            defaults to the outlet of a BlowdownMixer facility (if any).
        operating_hours : float, optional
            Number of operating hours in a year. This parameter is used to
            compute convinience properties such as utility cost and material cost
            on a per year basis.
        lang_factor : float, optional
            Lang factor for getting fixed capital investment from
            total purchase cost. If no lang factor, installed equipment costs are
            estimated using bare module factors.

        """
        if units is None:
            units = ()
        elif feeds is None:
            isa = isinstance
            Facility = bst.Facility
            feeds = bst.utils.feeds_from_units([i for i in units if not isa(i, Facility)])
            bst.utils.sort_feeds_big_to_small(feeds)
        if feeds:
            feedstock, *feeds = feeds
            facilities = facilities_from_units(units) if units else ()
            if not ends:
                ends = bst.utils.products_from_units(units) + [i.get_stream() for i in disjunctions]
            system = cls.from_feedstock(
                ID, feedstock, feeds, facilities, ends,
                facility_recycle or find_blowdown_recycle(facilities),
                operating_hours=operating_hours, lang_factor=lang_factor,
            )
        else:
            system = cls(ID, (), operating_hours=operating_hours)
        return system

    @classmethod
    def from_network(cls, ID, network, facilities=(), facility_recycle=None,
                     operating_hours=None, lang_factor=None):
        """
        Create a System object from a network.

        Parameters
        ----------
        ID : str
            Name of system.
        network : Network
            Network that defines the simulation path.
        facilities : Iterable[Facility]
            Offsite facilities that are simulated only after
            completing the path simulation.
        facility_recycle : [:class:`~thermosteam.Stream`], optional
            Recycle stream between facilities and system path.
        operating_hours : float, optional
            Number of operating hours in a year. This parameter is used to
            compute convinience properties such as utility cost and material cost
            on a per year basis.
        lang_factor : float, optional
            Lang factor for getting fixed capital investment from
            total purchase cost. If no lang factor, installed equipment costs are
            estimated using bare module factors.

        """
        facilities = Facility.ordered_facilities(facilities)
        isa = isinstance
        path = [(cls.from_network('', i) if isa(i, Network) else i)
                for i in network.path]
        self = cls.__new__(cls)
        self.recycle = network.recycle
        self._set_path(path)
        self._specification = None
        self._load_flowsheet()
        self._reset_errors()
        self._set_facilities(facilities)
        self._set_facility_recycle(facility_recycle)
        self._register(ID)
        self._load_defaults()
        self._save_configuration()
        self._load_stream_links()
        self.operating_hours = operating_hours
        self.lang_factor = lang_factor
        self._init_dynamic()
        return self

    def __init__(self, ID, path=(), recycle=None, facilities=(),
                 facility_recycle=None, N_runs=None, operating_hours=None,
                 lang_factor=None):
        self.recycle = recycle
        self.N_runs = N_runs
        self._set_path(path)
        self._specification = None
        self._load_flowsheet()
        self._reset_errors()
        self._set_facilities(facilities)
        self._set_facility_recycle(facility_recycle)
        self._register(ID)
        self._load_defaults()
        self._save_configuration()
        self._load_stream_links()
        self.operating_hours = operating_hours
        self.lang_factor = lang_factor
        self._init_dynamic()

    def __enter__(self):
        if self._path or self._recycle or self._facilities:
            raise RuntimeError("only empty systems can enter `with` statement")
        del self._unit_path, self._units, 
        unit_registry = self.flowsheet.unit
        self._irrelevant_units = set(unit_registry)
        unit_registry._open_dump(self)
        return self

    def __exit__(self, type, exception, traceback):
        irrelevant_units = self._irrelevant_units
        ID = self._ID
        del self._irrelevant_units
        unit_registry = self.flowsheet.unit
        dump = unit_registry._close_dump(self)
        if exception: raise exception
        if self._path or self._recycle or self._facilities:
            raise RuntimeError('system cannot be modified before exiting `with` statement')
        else:
            units = [i for i in dump if i not in irrelevant_units]
            system = self.from_units(None, units)
            self.ID = ID
            self.copy_like(system)

    def _save_configuration(self):
        self._connections = [i.get_connection() for i in bst.utils.streams_from_units(self.unit_path)]

    @ignore_docking_warnings
    def _load_configuration(self):
        for i in self._connections:
            if i.source:
                i.source.outs[i.source_index] = i.stream
            if i.sink:
                i.sink.ins[i.sink_index] = i.stream

    @ignore_docking_warnings
    def _interface_property_packages(self):
        path = self._path
        Stream = bst.Stream
        Interface = (bst.Junction, bst.Mixer, bst.MixTank)
        isa = isinstance
        new_path = []
        for obj in path:
            new_path.append(obj)
            outs = obj.outs
            for s in outs:
                source = s._source
                sink = s._sink
                if not sink or isa(sink, Interface): continue
                if sink.chemicals is not source.chemicals:
                    chemicals = s.chemicals
                    source_index = source._outs.index(s)
                    sink_index = sink._ins.index(s)
                    if chemicals is sink.chemicals:
                        s_sink = s
                        s_source = Stream(thermo=source.thermo)
                        s_source.copy_like(s)
                    else:
                        s_sink = Stream(thermo=sink.thermo)
                        s_sink.copy_like(s)
                        if chemicals is source.chemicals:
                            s_source = s
                        else:
                            s_source = Stream(thermo=source.thermo)
                            s_source.copy_like(s)
                    junction = bst.Junction(upstream=s_source, downstream=s_sink)
                    new_path.append(junction)
                    source._outs[source_index] = s_source
                    sink._ins[sink_index] = s_sink
        for obj in path:
            if isa(obj, System): obj._interface_property_packages()
        self._path = tuple(new_path)
        self._save_configuration()

    def _reduced_thermo_data(self, required_chemicals, unit_thermo, mixer_thermo, thermo_cache):
        isa = isinstance
        mixers = [i for i in self.units if isa(i, (bst.Mixer, bst.MixTank))]
        past_upstream_units = set()
        for mixer in mixers:
            if mixer in past_upstream_units: continue
            upstream_units = mixer.get_upstream_units()
            upstream_units.difference_update(past_upstream_units)
            available_chemicals = set(required_chemicals)
            for unit in upstream_units:
                if isa(unit, bst.Junction): continue
                available_chemicals.update(unit.get_available_chemicals())
            for unit in upstream_units:
                if isa(unit, bst.Junction): continue
                chemicals = [i for i in unit.chemicals if i in available_chemicals]
                if unit in unit_thermo:
                    other_thermo = unit_thermo[unit]
                    for i in other_thermo.chemicals:
                        if i not in chemicals: chemicals.append(i)
                IDs = tuple([i.ID for i in chemicals])
                if IDs in thermo_cache:
                    unit_thermo[unit] = thermo_cache[IDs]
                else:
                    unit_thermo[unit] = thermo_cache[IDs] = unit.thermo.subset(chemicals)
            past_upstream_units.update(upstream_units)
        for mixer in mixers:
            outlet = mixer.outs[0]
            sink = outlet.sink
            if sink:
                chemicals = sink.thermo.chemicals
            else:
                chemicals = outlet.available_chemicals
            if mixer in mixer_thermo:
                other_thermo = mixer_thermo[mixer]
                new_chemicals = []
                for i in other_thermo.chemicals:
                    if i not in chemicals: new_chemicals.append(i)
                if new_chemicals:
                    chemicals = list(chemicals) + new_chemicals
            IDs = tuple([i.ID for i in chemicals])
            if IDs in thermo_cache:
                mixer_thermo[mixer] = thermo_cache[IDs]
            else:
                mixer_thermo[mixer] = thermo_cache[IDs] = unit.thermo.subset(chemicals)

    def _delete_path_cache(self):
        for i in ('_units', '_unit_path', '_streams'):
            if hasattr(self, i): delattr(self, i)
        for i in self.subsystems: i._delete_path_cache()

    def reduce_chemicals(self, required_chemicals=()):
        self._delete_path_cache()
        unit_thermo = {}
        mixer_thermo = {}
        thermo_cache = {}
        self._reduced_thermo_data(required_chemicals, unit_thermo, mixer_thermo, thermo_cache)
        for unit, thermo in unit_thermo.items(): unit._reset_thermo(thermo)
        for mixer, thermo in mixer_thermo.items():
            for i in mixer._ins:
                if i._source: i._reset_thermo(unit_thermo[i._source])
            thermo = mixer_thermo[mixer]
            mixer._load_thermo(thermo)
            mixer._outs[0]._reset_thermo(thermo)
        self._interface_property_packages()

    def copy(self, ID=None):
        """Copy system."""
        new = System(ID)
        new.copy_like(self)
        return new

    def copy_like(self, other):
        """Copy path, facilities and recycle from other system."""
        self._path = other._path
        self._subsystems = other._subsystems
        self._facilities = other._facilities
        self._facility_loop = other._facility_loop
        self._recycle = other._recycle
        self._connections = other._connections

    def set_tolerance(self, mol=None, rmol=None, T=None, rT=None, subsystems=False, maxiter=None):
        """
        Set the convergence tolerance of the system.

        Parameters
        ----------
        mol : float, optional
            Molar tolerance.
        rmol : float, optional
            Relative molar tolerance.
        T : float, optional
            Temperature tolerance.
        rT : float, optional
            Relative temperature tolerance.
        subsystems : bool, optional
            Whether to also set tolerance of subsystems as well.
        maxiter : int, optional
            Maximum number if iterations.

        """
        if mol: self.molar_tolerance = float(mol)
        if rmol: self.relative_molar_tolerance = float(rmol)
        if T: self.temperature_tolerance = float(T)
        if rT: self.temperature_tolerance = float(rT)
        if maxiter: self.maxiter = int(maxiter)
        if subsystems:
            for i in self.subsystems: i.set_tolerance(mol, rmol, T, rT, subsystems, maxiter)

    ins = MockSystem.ins
    outs = MockSystem.outs
    load_inlet_ports = MockSystem.load_inlet_ports
    load_outlet_ports = MockSystem.load_outlet_ports
    _load_flowsheet  = MockSystem._load_flowsheet

    def _load_stream_links(self):
        for u in self.units: u._load_stream_links()

    def _load_defaults(self):
        #: [int] Maximum number of iterations.
        self.maxiter = self.default_maxiter

        #: [float] Molar tolerance (kmol/hr)
        self.molar_tolerance = self.default_molar_tolerance

        #: [float] Relative molar tolerance
        self.relative_molar_tolerance = self.default_relative_molar_tolerance

        #: [float] Temperature tolerance (K)
        self.temperature_tolerance = self.default_temperature_tolerance

        #: [float] Relative temperature tolerance
        self.relative_temperature_tolerance = self.default_relative_temperature_tolerance

        #: [str] Converge method
        self.converge_method = self.default_converge_method

        self.use_stabilized_convergence_algorithm = self.default_stabilized_convergence

    @property
    def TEA(self):
        """TEA object linked to the system."""
        return getattr(self, '_TEA', None)

    @property
    def LCA(self):
        """TEA object linked to the system."""
        return getattr(self, '_LCA', None)

    @property
    def specification(self):
        """Process specification."""
        return self._specification
    @specification.setter
    def specification(self, specification):
        if specification:
            if callable(specification):
                self._specification = specification
            else:
                raise AttributeError(
                    "specification must be callable or None; "
                   f"not a '{type(specification).__name__}'"
                )
        else:
            self._specification = None

    @property
    def use_stabilized_convergence_algorithm(self):
        """[bool] Whether to use a stablized convergence algorithm that implements 
        an inner loop with mass and energy balance approximations when applicable."""
        return self._stabilized
    @use_stabilized_convergence_algorithm.setter
    def use_stabilized_convergence_algorithm(self, stabilized):
        if stabilized and not self._recycle:
            for i in self.subsystems: i.use_stabilized_convergence_algorithm = True
        else:
            for i in self.subsystems: i.use_stabilized_convergence_algorithm = False
        self._stabilized = stabilized

    save_report = save_report

    def _extend_recycles(self, recycles):
        isa = isinstance
        recycle = self._recycle
        if recycle:
            if isa(recycle, Stream):
                recycles.append(recycle)
            elif isa(recycle, Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_error(recycle)
        for i in self._path:
            if isa(i, System): i._extend_recycles(recycles)

    def get_all_recycles(self):
        recycles = []
        self._extend_recycles(recycles)
        return recycles

    def _extend_flattend_path_and_recycles(self, path, recycles, stacklevel):
        isa = isinstance
        recycle = self._recycle
        stacklevel += 1
        if recycle:
            if isa(recycle, Stream):
                recycles.append(recycle)
            elif isa(recycle, Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_error(recycle)
        for i in self._path:
            if isa(i, System):
                if i.facilities:
                    warning = RuntimeWarning('subsystem with facilities could not be flattened')
                    warn(warning, stacklevel=stacklevel)
                    path.append(i)
                elif i.specification:
                    warning = RuntimeWarning('subsystem with specification could not be flattened')
                    warn(warning, stacklevel=stacklevel)
                    path.append(i)
                else:
                    i._extend_flattend_path_and_recycles(path, recycles, stacklevel)
            else:
                path.append(i)

    def prioritize_unit(self, unit):
        """
        Prioritize unit operation to run first within it's recycle system,
        if there is one.

        Parameters
        ----------
        unit : Unit
            Unit operation to prioritize.

        Raises
        ------
        ValueError
            When unit is not in the system.
        RuntimeError
            When prioritization algorithm fails. This should never happen.

        Examples
        --------
        Create a simple recycle loop and prioritize a different unit operation:

        >>> from biosteam import main_flowsheet as f, Stream, settings, Mixer, Splitter
        >>> f.set_flowsheet('simple_recycle_loop')
        >>> settings.set_thermo(['Water'], cache=True)
        >>> feedstock = Stream('feedstock', Water=1000)
        >>> water = Stream('water', Water=10)
        >>> recycle = Stream('recycle')
        >>> product = Stream('product')
        >>> M1 = Mixer('M1', [feedstock, water, recycle])
        >>> S1 = Splitter('S1', M1-0, [product, recycle], split=0.5)
        >>> recycle_loop_sys = f.create_system('recycle_loop_sys')
        >>> recycle_loop_sys.print()
        System('recycle_loop_sys',
            [M1,
             S1],
            recycle=S1-1)
        >>> recycle_loop_sys.prioritize_unit(S1)
        >>> recycle_loop_sys.print()
        System('recycle_loop_sys',
            [S1,
             M1],
            recycle=S1-1)

        """
        isa = isinstance
        if unit not in self.unit_path:
            raise ValueError(f'unit {repr(unit)} not in system')
        path = self._path
        if (self._recycle or self.N_runs):
            for index, other in enumerate(path):
                if unit is other:
                    self._path = path[index:] + path[:index]
                    return
                elif isa(other, System) and unit in other.unit_path:
                    other.prioritize_unit(unit)
                    return
            raise RuntimeError('problem in system algorithm')


    def split(self, stream, ID_upstream=None, ID_downstream=None):
        """
        Split system in two; upstream and downstream.

        Parameters
        ----------
        stream : Iterable[:class:~thermosteam.Stream], optional
            Stream where unit group will be split.
        ID_upstream : str, optional
            ID of upstream system.
        ID_downstream : str, optional
            ID of downstream system.

        Examples
        --------
        >>> from biorefineries.cornstover import cornstover_sys, M201
        >>> from biosteam import default
        >>> upstream_sys, downstream_sys = cornstover_sys.split(M201-0)
        >>> upstream_group = upstream_sys.to_unit_group()
        >>> upstream_group.show()
        UnitGroup: Unnamed
         units: U101, H2SO4_storage, T201, M201
        >>> downstream_group = downstream_sys.to_unit_group()
        >>> for i in upstream_group: assert i not in downstream_group.units
        >>> assert set(upstream_group.units + downstream_group.units) == set(cornstover_sys.units)
        >>> default() # Reset to biosteam defaults

        """
        if self._recycle: raise RuntimeError('cannot split system with recycle')
        path = self._path
        streams = self.streams
        surface_units = {i for i in path if isinstance(i, Unit)}
        if stream.source in surface_units:
            index = path.index(stream.source) + 1
        elif stream.sink in surface_units:
            index = path.index(stream.sink)
        elif stream not in streams:
            raise ValueError('stream not in system')
        else:
            raise ValueError('stream cannot reside within a subsystem')
        return (System(ID_upstream, path[:index], None),
                System(ID_downstream, path[index:], None, self._facilities))

    def flatten(self):
        """Flatten system by removing subsystems."""
        recycles = []
        path = []
        self._extend_flattend_path_and_recycles(path, recycles, stacklevel=2)
        self._path = tuple(path)
        self._recycle = tuple(recycles)
        N_recycles = len(recycles)
        self.molar_tolerance *= N_recycles
        self.temperature_tolerance *= N_recycles

    def to_unit_group(self, name=None):
        """Return a UnitGroup object of all units within the system."""
        return bst.UnitGroup(name, self.units)

    def _set_path(self, path):
        #: tuple[Unit, function and/or System] A path that is run element
        #: by element until the recycle converges.
        self._path = path = tuple(path)

    def _set_facilities(self, facilities):
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the path simulation.
        self._facilities = tuple(facilities)
        self._load_facilities()

    def _load_facilities(self):
        isa = isinstance
        units = self.units.copy()
        for i in self._facilities:
            if isa(i, Facility):
                if i._system: continue
                i._system = self
                i._other_units = other_units = units.copy()
                other_units.remove(i)


    def _set_facility_recycle(self, recycle):
        if recycle:
            sys = self._downstream_system(recycle._sink)
            sys.recycle = recycle
            sys.__class__ = FacilityLoop
            #: [FacilityLoop] Recycle loop for converging facilities
            self._facility_loop = sys
        else:
            self._facility_loop = None

    # Forward pipping
    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__

    @property
    def subsystems(self):
        """list[System] All subsystems in the system."""
        try:
            return self._subsystems
        except:
            self._subsystems = [i for i in self._path if isinstance(i, System)]
            return self._subsystems

    @property
    def units(self):
        """[list] All unit operations as ordered in the path without repetitions."""
        try:
            return self._units
        except:
            self._units = units = []
            past_units = set()
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit):
                    if i in past_units: continue
                    units.append(i)
                    past_units.add(i)
                elif isa(i, System):
                    sys_units = i.units
                    units.extend([i for i in sys_units if i not in past_units])
                    past_units.update(sys_units)
            return units

    @property
    def unit_path(self):
        """[list] Unit operations as ordered in the path (some units may be repeated)."""
        try:
            return self._unit_path
        except:
            self._unit_path = units = []
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit):
                    units.append(i)
                elif isa(i, System):
                    units.extend(i.unit_path)
            return units

    @property
    def cost_units(self):
        """[set] All unit operations with costs."""
        try:
            return self._cost_units
        except:
            self._cost_units = units = set()
            isa = isinstance
            for i in self._path + self._facilities:
                if isa(i, Unit) and (i._design or i._cost):
                    units.add(i)
                elif isa(i, System):
                    units.update(i.cost_units)
            return units
        
    @property
    def streams(self):
        """set[:class:`~thermosteam.Stream`] All streams within the system."""
        try:
            return self._streams
        except:
            self._streams = streams = bst.utils.streams_from_units(self.unit_path)
            bst.utils.filter_out_missing_streams(streams)
            return streams
    @property
    def feeds(self):
        """set[:class:`~thermosteam.Stream`] All feeds to the system."""
        try:
            return self._feeds
        except:
            self._feeds = feeds = bst.utils.feeds(self.streams)
            return feeds
    @property
    def products(self):
        """set[:class:`~thermosteam.Stream`] All products of the system."""
        try:
            return self._products
        except:
            self._products = products = bst.utils.products(self.streams)
            return products

    @property
    def facilities(self):
        """tuple[Facility] All system facilities."""
        return self._facilities

    @property
    def recycle(self):
        """
        :class:`~thermosteam.Stream` or Iterable[:class:`~thermosteam.Stream`]
        A tear stream for the recycle loop.
        """
        return self._recycle
    @recycle.setter
    def recycle(self, recycle):
        isa = isinstance
        self._N_runs = None
        if recycle is None:
            self._recycle = recycle
        elif isa(recycle, Stream):
            self._recycle = recycle
        elif isa(recycle, Iterable):
            recycle = set(recycle)
            for i in recycle:
                if not isa(i, Stream):
                    raise ValueError("recycle streams must be Stream objects; "
                                     f"not {type(i).__name__}")
            self._recycle = recycle
        else:
            raise_recycle_type_error(recycle)

    @property
    def N_runs(self):
        """Number of times to run the path."""
        return self._N_runs
    @N_runs.setter
    def N_runs(self, N_runs):
        if N_runs: self._recycle = None
        self._N_runs = N_runs

    @property
    def path(self):
        """
        tuple[Unit, function and/or System] A path that is run element by
        element until the recycle(s) converges (if any).

        """
        return self._path

    @property
    def converge_method(self):
        """Iterative convergence method ('wegstein', 'aitken', or 'fixedpoint')."""
        return self._converge_method.__name__[1:]
    @converge_method.setter
    def converge_method(self, method):
        method = method.lower().replace('-', '').replace(' ', '')
        try:
            self._converge_method = getattr(self, '_' + method)
        except:
            raise ValueError("only 'wegstein', 'aitken', and 'fixedpoint' "
                            f"methods are valid, not '{method}'")

    @property
    def isdynamic(self):
        '''Whether the system contains any dynamic Unit.'''
        isdynamic = [(unit._isdynamic if hasattr(unit, '_isdynamic') else False)
                     for unit in self.units]
        return any(isdynamic)

    def _downstream_path(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.unit_path: return []
        elif self._recycle: return self._path
        isa = isinstance
        for index, i in enumerate(self._path):
            if unit is i:
                return self._path[index:]
            elif (isa(i, System) and unit in i.unit_path):
                return i._downstream_path(unit) + self._path[index+1:]
        return []

    def _downstream_facilities(self, unit):
        """Return a list of facilities composed of the `unit` and
        everything downstream."""
        isa = isinstance
        for index, i in enumerate(self._facilities):
            if unit is i or (isa(i, System) and unit in i.unit_path):
                return self._facilities[index:]
        return []

    def _downstream_system(self, unit):
        """Return a system with a path composed of the `unit` and
        everything downstream (facilities included)."""
        if self._recycle or unit is self._path[0]: return self
        path = self._downstream_path(unit)
        if path:
            facilities = self._facilities
        else:
            facilities = self._downstream_facilities(unit)
            if not facilities:
                raise RuntimeError(f'{unit} not found in system')
        system = System(None, path,
                        facilities=facilities)
        system._ID = f'{type(unit).__name__}-{unit} and downstream'
        return system

    def _minimal_digraph(self, graph_attrs):
        """Return digraph of the path as a box."""
        return minimal_digraph(self.ID, self.units, self.streams, **graph_attrs)

    def _surface_digraph(self, graph_attrs):
        return surface_digraph(self._path, **graph_attrs)

    def _thorough_digraph(self, graph_attrs):
        return digraph_from_units_and_streams(self.unit_path,
                                              self.streams,
                                              **graph_attrs)

    def _cluster_digraph(self, graph_attrs):
        return digraph_from_system(self, **graph_attrs)

    def diagram(self, kind=None, file=None, format=None, display=True,
                number=None, profile=None, label=None, **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of
        the system.

        Parameters
        ----------
        kind : int or string, optional
            * 0 or 'cluster': Display all units clustered by system.
            * 1 or 'thorough': Display every unit within the path.
            * 2 or 'surface': Display only elements listed in the path.
            * 3 or 'minimal': Display a single box representing all units.
        file=None : str, display in console by default
            File name to save diagram.
        format='png' : str
            File format (e.g. "png", "svg").
        display : bool, optional
            Whether to display diagram in console or to return the graphviz
            object.
        number : bool, optional
            Whether to number unit operations according to their
            order in the system path.
        profile : bool, optional
            Whether to clock the simulation time of unit operations.
        label : bool, optional
            Whether to label the ID of streams with sources and sinks.

        """
        self._load_configuration()
        if not kind: kind = 0
        graph_attrs['format'] = format or 'png'
        original = (bst.LABEL_PATH_NUMBER_IN_DIAGRAMS,
                    bst.LABEL_PROCESS_STREAMS_IN_DIAGRAMS,
                    bst.PROFILE_UNITS_IN_DIAGRAMS)
        if number is not None: bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = number
        if label is not None: bst.LABEL_PROCESS_STREAMS_IN_DIAGRAMS = label
        if profile is not None: bst.PROFILE_UNITS_IN_DIAGRAMS = profile
        try:
            if kind == 0 or kind == 'cluster':
                f = self._cluster_digraph(graph_attrs)
            elif kind == 1 or kind == 'thorough':
                f = self._thorough_digraph(graph_attrs)
            elif kind == 2 or kind == 'surface':
                f = self._surface_digraph(graph_attrs)
            elif kind == 3 or kind == 'minimal':
                f = self._minimal_digraph(graph_attrs)
            else:
                raise ValueError("kind must be one of the following: "
                                 "0 or 'cluster', 1 or 'thorough', 2 or 'surface', "
                                 "3 or 'minimal'")
            if display or file:
                finalize_digraph(f, file, format)
            else:
                return f
        finally:
            (bst.LABEL_PATH_NUMBER_IN_DIAGRAMS,
             bst.LABEL_PROCESS_STREAMS_IN_DIAGRAMS,
             bst.PROFILE_UNITS_IN_DIAGRAMS) = original

    # Methods for running one iteration of a loop
    def _iter_run(self, mol):
        """
        Run the system at specified recycle molar flow rate.

        Parameters
        ----------
        mol : numpy.ndarray
              Recycle molar flow rates.

        Returns
        -------
        mol_new : numpy.ndarray
            New recycle molar flow rates.
        not_converged : bool
            True if recycle has not converged.

        """
        mol[mol < 0.] = 0.
        self._set_recycle_data(mol)
        T = self._get_recycle_temperatures()
        self._run()
        mol_new = self._get_recycle_data()
        T_new = self._get_recycle_temperatures()
        mol_errors = np.abs(mol - mol_new)
        positive_index = mol_errors > 1e-16
        mol_errors = mol_errors[positive_index]
        if mol_errors.size == 0:
            self._mol_error = mol_error = 0.
            self._rmol_error = rmol_error = 0.
        else:
            self._mol_error = mol_error = mol_errors.max()
            if mol_error > 1e-12:
                self._rmol_error = rmol_error = (mol_errors / np.maximum.reduce([np.abs(mol[positive_index]), np.abs(mol_new[positive_index])])).max()
            else:
                self._rmol_error = rmol_error = 0.
        T_errors = np.abs(T - T_new)
        self._T_error = T_error = T_errors.max()
        self._rT_error = rT_error = (T_errors / T).max()
        self._iter += 1
        not_converged = not (
            (mol_error < self.molar_tolerance
             or rmol_error < self.relative_molar_tolerance)
            and
            (T_error < self.temperature_tolerance
             or rT_error < self.relative_temperature_tolerance)
        )
        if not_converged and self._iter >= self.maxiter:
            if self.strict_convergence: raise RuntimeError(f'{repr(self)} could not converge' + self._error_info())
            else: not_converged = False
        return mol_new, not_converged

    def _get_recycle_data(self):
        recycle = self._recycle
        if isinstance(recycle, Stream):
            return recycle._imol.data.copy()
        elif isinstance(recycle, Iterable):
            return np.vstack([i._imol.data for i in recycle])
        else:
            raise RuntimeError('no recycle available')

    def _set_recycle_data(self, data):
        recycle = self._recycle
        isa = isinstance
        if isa(recycle, Stream):
            try:
                recycle._imol._data[:] = data
            except IndexError as e:
                if data.shape[0] != 1:
                    raise IndexError(f'expected 1 row; got {data.shape[0]} rows instead')
                else:
                    raise e from None
        elif isa(recycle, Iterable):
            length = len
            N_rows = data.shape[0]
            M_rows = sum([length(i) if isa(i, MultiStream) else 1 for i in recycle])
            if M_rows != N_rows:
                raise IndexError(f'expected {M_rows} rows; got {N_rows} rows instead')
            index = 0
            for i in recycle:
                if isa(i, MultiStream):
                    next_index = index + length(i)
                    i._imol._data[:] = data[index:next_index, :]
                    index = next_index
                else:
                    i._imol._data[:] = data[index, :]
                    index += 1
        else:
            raise RuntimeError('no recycle available')

    def _get_recycle_temperatures(self):
        recycle = self._recycle
        if isinstance(recycle, Stream):
            T = self._recycle.T
        elif isinstance(recycle, Iterable):
            T = [i.T for i in recycle]
        else:
            raise RuntimeError('no recycle available')
        return np.array(T, float)

    def _setup(self):
        """Setup each element of the system."""
        self._load_facilities()
        self._load_configuration()
        for i in self.units: i._setup()

    def _run(self):
        """Rigorously run each element in the path."""
        isa = isinstance
        converge = converge_system_in_path
        run = try_method_with_object_stamp
        for i in self._path:
            if isa(i, Unit): run(i, i.run)
            elif isa(i, System): converge(i)
            else: i() # Assume it's a function

    # Methods for convering the recycle stream
    def _fixedpoint(self):
        """Converge system recycle iteratively using fixed-point iteration."""
        self._solve(flx.conditional_fixed_point)

    def _wegstein(self):
        """Converge the system recycle iteratively using wegstein's method."""
        self._solve(flx.conditional_wegstein)

    def _aitken(self):
        """Converge the system recycle iteratively using Aitken's method."""
        self._solve(flx.conditional_aitken)

    def _solve(self, solver):
        """Solve the system recycle iteratively using given solver."""
        self._reset_iter()
        f = iter_run = self._iter_run
        if self._stabilized:
            special_units = [i for i in self.units if hasattr(i, '_steady_run')]
            if special_units:
                def f(mol):
                    self._set_recycle_data(mol)
                    for unit in special_units: unit._run = unit._steady_run
                    try:
                        solver(iter_run, self._get_recycle_data())
                    finally:
                        for unit in special_units: del unit._run
                    return iter_run(self._get_recycle_data())
        try:
            solver(f, self._get_recycle_data())
        except IndexError as error:
            try:
                solver(f, self._get_recycle_data())
            except:
                raise error

    def _converge(self):
        if self._N_runs:
            for i in range(self.N_runs): self._run()
        elif self._recycle:
            self._converge_method()
        else:
            self._run()

    def _summary(self):
        simulated_units = set()
        isa = isinstance
        Unit = bst.Unit
        for i in self._path:
            if isa(i, Unit):
                if i in simulated_units: continue
                simulated_units.add(i)
            try_method_with_object_stamp(i, i._summary)
        simulate_unit = simulate_unit_in_path
        for i in self._facilities:
            if isa(i, Unit): simulate_unit(i)
            elif isa(i, System):
                converge_system_in_path(i)
                i._summary()
            else: i() # Assume it is a function
        for i in self._facilities:
            if isa(i, (bst.BoilerTurbogenerator, bst.Boiler)): simulate_unit(i)

    def _reset_iter(self):
        self._iter = 0
        for system in self.subsystems: system._reset_iter()

    def _reset_errors(self):
        #: Molar flow rate error (kmol/hr)
        self._mol_error = 0

        #: Relative molar flow rate error
        self._rmol_error = 0

        #: Temperature error (K)
        self._T_error = 0

        #: Relative temperature error
        self._rT_error = 0

        #: Number of iterations
        self._iter = 0

    def empty_outlet_streams(self):
        """Reset all outlet streams to zero flow."""
        self._reset_errors()
        units = self.units
        streams = bst.utils.streams_from_units(units)
        bst.utils.filter_out_missing_streams(streams)
        streams_by_data = {}
        for i in streams:
            data = i.imol.data
            data_id = id(data)
            if data_id in streams_by_data:
                streams_by_data[data_id].append(i)
            else:
                streams_by_data[data_id] = [i]
        for streams in streams_by_data.values():
            if all([i.source in units for i in streams]):
                streams[0].empty()

    def empty_recycles(self):
        """Reset all recycle streams to zero flow."""
        self._reset_errors()
        recycle = self._recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycle.empty()
            elif isinstance(recycle, Iterable):
                for i in recycle: i.empty()
            else:
                raise_recycle_type_error(recycle)
        for system in self.subsystems:
            system.empty_recycles()

    def _init_dynamic(self):
        '''Initialize attributes related to dynamic simulation.'''
        self._state = None

    def reset_cache(self):
        """Reset cache of all unit operations."""
        for unit in self.units: unit.reset_cache()
        self._init_dynamic()

    def _state_attr2arr(self):
        arr = np.array([])
        idxer = {}
        for unit in self.units:
            start = len(arr)
            arr = np.append(arr, unit._state)
            stop = len(arr)
            idxer[unit._ID] = (start, stop)
        return arr, idxer

    def _dstate_attr2arr(self, arr, idx):
        dy = np.zeros_like(arr)
        for unit in self.units:
            start, stop = idx[unit._ID]
            dy[start: stop] = unit._dstate
        return dy

    def _state_arr2attr(self, arr, idx):
        for unit in self.units:
            start, stop = idx[unit._ID]
            unit._update_state(arr[start: stop])

    def _load_state(self):
        '''Returns the initial state (a 1d-array) of the system for dynamic simulation.'''
        if self._state is None:
            for ws in self.feeds:
                if ws._state is None: ws._init_state()
            n_rotate = 0
            for u in self.units:
                if not u.isdynamic: n_rotate += 1
                else: break
            units = self.units[n_rotate:] + self.units[:n_rotate]
            for inf in units[0].ins:
                inf._init_state()
            for unit in units: 
                try:
                    if unit._state is None: unit._init_state()   
                except AttributeError:
                    raise AttributeError(f'{unit.ID} does not have `_init_state`, check if it is a dynamic unit.')
                unit._update_state(unit._state)
                unit._update_dstate()  
            y, idx = self._state_attr2arr()
            self._state = {'time': 0,
                           'state': y,
                           'indexer': idx,
                           'state_over_time': None,
                           'n_rotate': n_rotate}
        else:
            y = self._state['state']
            idx = self._state['indexer']
            n_rotate = self._state['n_rotate']
            units = self.units[n_rotate:] + self.units[:n_rotate]
            for unit in units: unit._update_dstate()  
        return y, idx, n_rotate

    def _ODE(self, idx, n_rotate):
        '''System-wide ODEs.'''
        units = self.units[n_rotate:] + self.units[:n_rotate]
        def dydt(t, y):
            self._state_arr2attr(y, idx)
            for unit in units:
                QC_ins = unit._collect_ins_state()
                dQC_ins = unit._collect_ins_dstate()
                QC = unit._state
                unit.ODE(t, QC_ins, QC, dQC_ins)   
            return self._dstate_attr2arr(y, idx)
        return dydt

    def _write_state(self, t, y):
        '''record the states at a certain time point of the dynamic simulation and define wastestreams '''
        self._state['time'] = t
        self._state['state'] = y
        idx = self._state['indexer']
        self._state_arr2attr(y, idx)
        for ws in self.streams - set(self.feeds):
            ws._state2flows()

    def clear_state(self):
        self._state = None
        for u in self.units:
            u._state = None
        for ws in self.streams:
            ws._state = None

    def _state2df(self, path=''):
        header = [('-', 't [d]')] + sum([[(m, n) for m, n in zip([u.ID]*len(u._state_header), u._state_header)] for u in self.units], [])
        import pandas as pd
        header = pd.MultiIndex.from_tuples(header, names=['unit', 'variable'])
        data = self._state['state_over_time']
        df = pd.DataFrame(data, columns=header)
        if not path:
            return df
        if path.endswith(('.xlsx', '.xls')):
            df.to_excel(path)
        elif path.endswith('.csv'):
            df.to_csv(path)
        elif path.endswith('.tsv'):
            df.to_csv(path, sep='\t')
        else:
            ext = path.split('.')[-1]
            raise ValueError('Only support file extensions of '
                             '".xlsx", ".xls", ".csv", and ".tsv", '
                             f'not .{ext}.')

    def simulate(self, start_from_cached_state=True,
                 solver='solve_ivp', export_state_to='', print_msg=False,
                 **kwargs):
        """
        Converge the path and simulate all units.
        
        Parameters
        ----------
        start_from_cached_state: bool
            Whether to start from the cached state.
        solver : str
            Which ``scipy`` function to use, either "solve_ivp" or "odeint".
        export_state_to: str
            If provided with a path, will save the simulated states over time to the given path,
            supported extensions are '.xlsx', '.xls', 'csv', and 'tsv'.
        kwargs : dict
            Other keyword arguments that will be passed to ``solve_ivp``
            or ``odeint``. Must contain t_span for ``solve_ivp`` or t for ``odeint``.
        
        See Also
        --------
        `scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_
        `scipy.integrate.odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_        
        """
        self._setup()
        self._converge()
        if self.isdynamic:
            if not start_from_cached_state: self.clear_state()
            y0, idx, nr = self._load_state()
            dydt = self._ODE(idx, nr)
            if solver == 'solve_ivp':
                sol = solve_ivp(fun=dydt, y0=y0, **kwargs)
                self._state['state_over_time'] = np.vstack((sol.t, sol.y)).T
                if sol.status == 0:
                    print('Simulation completed.')
                else:
                    print(sol.message)
                self._write_state(sol.t[-1], sol.y.T[-1])
            elif solver=='odeint':
                sol = odeint(func=dydt, y0=y0, printmessg=print_msg, tfirst=True, **kwargs)
                t_arr = kwargs['t']
                if sol.shape[0] < len(t_arr):
                    print('Simulation failed.')
                else:
                    print('Simulation completed.')
                self._state['state_over_time'] = np.hstack((t_arr.reshape((len(t_arr), 1)), sol))
                self._write_state(t_arr[-1], sol[-1])
            else:
                raise ValueError('`solver` can only be "solve_ivp" or "odeint", '
                                 f'not {solver}.')
            if export_state_to:
                self._state2df(export_state_to)
        self._summary()
        if self._facility_loop: self._facility_loop._converge()

    # User definitions
    
    def define_process_impact(self, key, name, basis, inventory, CF):
        """
        Define a process impact.

        Parameters
        ----------
        key : str
            Impact indicator key.
        name : str
            Name of process impact.
        basis : str
            Functional unit for the charaterization factor.
        inventory : callable
            Should return the annualized (not hourly) inventory flow rate 
            given no parameters. Units should be in basis / yr
        CF : float
            Characterization factor in impact / basis.

        """
        item = ProcessImpactItem(name, basis, inventory, CF)
        try:
            if key in self.process_impact_items:
                items = self.process_impact_items[key]
                name = item.name
                for i in items:
                    if i.name == name:
                        raise ValueError(
                            f"item with name '{name}' already defined"
                        )
                items.append(item)
            else:
                self.process_impact_items[key] = [item]
        except AttributeError:
            self.process_impact_items = {key: [item]}

    # Convenience methods

    @property
    def heat_utilities(self):
        """tuple[HeatUtility] the sum of all heat utilities in the system by agent."""
        return bst.HeatUtility.sum_by_agent(utils.get_heat_utilities(self.cost_units))

    @property
    def power_utility(self):
        """[PowerUtility] Sum of all power utilities in the system."""
        return bst.PowerUtility.sum(utils.get_power_utilities(self.cost_units))

    def get_inlet_utility_flows(self):
        """
        Return a dictionary with inlet stream utility flow rates, including
        natural gas and ash disposal but excluding heating, refrigeration, and
        electricity utilities.
        
        """
        dct = {}
        for unit in self.units:
            for name, flow in unit.get_inlet_utility_flows().items():
                if flow:
                    if name in dct:
                        dct[name] += flow
                    else:
                        dct[name] = flow
        return dct

    def get_outlet_utility_flows(self):
        """
        Return a dictionary with outlet stream utility flow rates, including
        natural gas and ash disposal but excluding heating, refrigeration, and
        electricity utilities.
        
        """
        dct = {}
        for unit in self.units:
            for name, flow in unit.get_outlet_utility_flows().items():
                if flow:
                    if name in dct:
                        dct[name] += flow
                    else:
                        dct[name] = flow
        return dct

    def get_inlet_flow(self, units, key=None):
        """
        Return total flow across all inlets per year.

        Parameters
        ----------
        units : str
            Material units of measure (e.g., 'kg', 'gal', 'kmol').
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned

        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> sys.get_inlet_flow('Mton') # Sum of all chemicals
        0.1584
        >>> sys.get_inlet_flow('Mton', 'Water') # Just water
        0.0792

        """
        units += '/hr'
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in bst.utils.feeds_from_units(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in bst.utils.feeds_from_units(self.units)])

    def get_outlet_flow(self, units, key=None):
        """
        Return total flow across all outlets per year.

        Parameters
        ----------
        units : str
            Material units of measure (e.g., 'kg', 'gal', 'kmol').
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned

        Examples
        --------
        >>> from biosteam import Stream, Mixer, Splitter, settings, main_flowsheet
        >>> settings.set_thermo(['Water', 'Ethanol'])
        >>> main_flowsheet.clear()
        >>> S1 = Splitter('S1', Stream(Ethanol=10, units='ton/hr'), split=0.1)
        >>> M1 = Mixer('M1', ins=[Stream(Water=10, units='ton/hr'), S1-0])
        >>> sys = main_flowsheet.create_system(operating_hours=330*24)
        >>> sys.simulate()
        >>> sys.get_outlet_flow('Mton') # Sum of all chemicals
        0.1584
        >>> sys.get_outlet_flow('Mton', 'Water') # Just water
        0.0792

        """
        units += '/hr'
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in bst.utils.products_from_units(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in bst.utils.products_from_units(self.units)])
    
    def get_mass_flow(self, stream):
        """Return the mass flow rate of a stream [kg/yr]."""
        return stream.F_mass * self.operating_hours
    
    def get_market_value(self, stream):
        """Return the market value of a stream [USD/yr]."""
        return stream.cost * self.operating_hours

    def get_property(self, stream, name, units=None):
        """Return the annualized property of a stream."""
        return stream.get_property(name, units) * self.operating_hours

    def _price2cost(self, stream):
        """Get factor to convert stream price to cost."""
        F_mass = stream.F_mass
        if not F_mass: warn(RuntimeWarning(f"stream '{stream}' is empty"))
        price2cost = F_mass * self.operating_hours
        if stream.sink and not stream.source:
            return - price2cost
        elif stream.source:
            return price2cost
        else:
            raise ValueError("stream must be either a feed or a product")

    def get_net_heat_utility_impact(self, agent, key, heat_utilities=None):
        if isinstance(agent, str): 
            ID = agent
            agent = bst.HeatUtility.get_agent(ID)
        elif isinstance(agent, bst.UtilityAgent):
            ID = agent.ID
        else:
            raise TypeError(
                 "agent must be a UtilityAgent object or a string; "
                f"not a '{type(agent).__name__}' object"
            )
        CF, units = bst.HeatUtility.get_CF(ID, key)
        if CF == 0.: return 0.
        if heat_utilities is None: heat_utilities = self.heat_utilities
        for hu in heat_utilities:
            if hu.agent and hu.agent.ID == ID:
                if units == 'kg':
                    return hu.flow * CF * agent.MW * self.operating_hours
                elif units == 'kmol':
                    return hu.flow * CF * self.operating_hours
                elif units == 'kJ':
                    return hu.duty * CF * self.operating_hours
                else:
                    raise RuntimeError("unknown error")
        return 0.
    
    def get_net_electricity_impact(self, key):
        try:
            return self.power_utility.get_impact(key) * self.operating_hours
        except KeyError:
            return 0.

    def get_net_utility_impact(self, key):
        agents = (*bst.HeatUtility.cooling_agents,
                  *bst.HeatUtility.heating_agents)
        heat_utilities = self.heat_utilities
        return sum([self.get_net_heat_utility_impact(i, key, heat_utilities) for i in agents]) * self.operating_hours + self.get_net_electricity_impact(key)
    
    def get_total_feeds_impact(self, key):
        """
        Return the total annual impact of all feeds given 
        the impact indicator key.
        
        """
        return sum([s.F_mass * s.characterization_factors[key] for s in self.feeds
                    if key in s.characterization_factors]) * self.operating_hours
    
    def get_total_products_impact(self, key):
        """
        Return the total annual impact of all products given 
        the impact indicator key.
        
        """
        return sum([s.F_mass * s.characterization_factors[key] for s in self.products
                    if key in s.characterization_factors]) * self.operating_hours
    
    def get_material_impact(self, stream, key):
        """
        Return the annual material impact given the stream and the 
        impact indicator key.
        
        """
        return stream.get_impact(key) * self.operating_hours
    
    def get_total_input_impact(self, key):
        """
        Return total annual impact of inputs given the impact indicator key.
        
        """
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        impact = self.get_total_feeds_impact(key) / self.operating_hours
        for hu in heat_utilities:
            if hu.flow > 0.: impact += hu.get_impact(key)
        if power_utility.rate > 0.:
            impact += power_utility.get_impact(key)
        return impact * self.operating_hours
    
    def get_process_impact(self, key):
        """
        Return the annual process impact given the impact indicator key.
        
        """
        try:
            process_impact_items = self.process_impact_items
        except:
            return 0.
        else:
            return sum([item.impact() for item in process_impact_items[key]])
    
    def get_net_impact(self, key):
        """
        Return net annual impact, including displaced impacts, given the impact indicator key.
        
        """
        return (
            self.get_total_feeds_impact(key)
            + self.get_process_impact(key)
            + self.get_net_utility_impact(key)
            - self.get_total_products_impact(key)
        )
    
    def get_property_allocated_impact(self, key, name, units):
        total_property = 0.
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        operating_hours = self.operating_hours
        if hasattr(bst.PowerUtility, name):
            if power_utility.rate < 0.:
                total_property += power_utility.get_property(name, units) * operating_hours
        if hasattr(bst.HeatUtility, name):
            for hu in heat_utilities:
                if hu.flow < 0.: total_property += hu.get_property(name, units) * operating_hours
        if hasattr(bst.Stream, name):
            for stream in self.products:
                total_property += self.get_property(stream, name, units)
        impact = self.get_total_feeds_impact(key)
        for hu in heat_utilities:
            if hu.flow > 0.: impact += hu.get_impact(key) * operating_hours
        if power_utility.rate > 0.:
            impact += power_utility.get_impact(key) * operating_hours
        impact += self.get_process_impact(key)
        return impact / total_property
    
    def get_property_allocation_factors(self, name, units):
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        operating_hours = self.operating_hours
        properties = {}
        if hasattr(bst.PowerUtility, name):
            if power_utility.rate < 0.:
                value = power_utility.get_property(name, units)
                if value: properties['Electricity'] = value * operating_hours
        if hasattr(bst.HeatUtility, name):
            for hu in heat_utilities:
                if hu.flow < 0.: 
                    value = hu.get_property(name, units)
                    if value: properties[_reformat(hu.agent.ID)] = value * operating_hours
        if hasattr(bst.Stream, name):
            for stream in self.products:
                value = self.get_property(stream, name, units)
                if value: properties[_reformat(stream.ID)] = value
        total_property = sum(properties.values())
        allocation_factors = {i: j / total_property for i, j in properties.items()}
        return allocation_factors
    
    def get_displacement_allocation_factors(self, main_product, key):
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        allocation_factors = {}
        isa = isinstance
        if isa(main_product, bst.Stream):
            CF_original = main_product.get_CF(key)
        elif main_product == 'Electricity':
            main_product = power_utility
            CF_original = main_product.get_CF(key, consumption=False)
        else:
            raise NotImplementedError(f"main product '{main_product}' is not yet an option for this method")
        items = [main_product]
        for stream in self.products:
            if key in stream.characterization_factors:
                items.append(stream)
        for hu in heat_utilities:
            if hu.flow > 0. and (hu.agent.ID, key) in hu.characterization_factors:
                items.append(hu)
        if power_utility.rate < 0. and key in power_utility.characterization_factors:
            items.append(power_utility)
        for item in items:
            total_input_impact = self.get_total_input_impact(key)
            net_impact = self.get_net_impact(key)
            if isa(item, bst.Stream):
                item_impact = self.get_material_impact(item, key)
            else:
                item_impact = -1 * item.get_impact(key) * self.operating_hours
            displaced_impact = net_impact + item_impact
            if isa(item, bst.Stream):
                allocation_factors[_reformat(item.ID)] = displaced_impact / total_input_impact
            elif isa(item, bst.HeatUtility):
                allocation_factors[_reformat(item.ID)] = displaced_impact / total_input_impact
            elif isa(item, bst.PowerUtility):
                allocation_factors['Electricity'] = displaced_impact / total_input_impact
            else:
                raise RuntimeError('unknown error')
            if item is main_product:
                if isa(main_product, bst.Stream):
                    main_product.set_CF(key, displaced_impact / self.get_mass_flow(main_product))
                elif main_product == 'Electricity':
                    main_product.set_CF(key, displaced_impact / (main_product.rate * self.operating_hours))
        main_product.set_CF(key, CF_original)
        total = sum(allocation_factors.values())
        return {i: j / total for i, j in allocation_factors.items()}
    
    @property
    def sales(self):
        """Annual sales revenue."""
        return sum([s.cost for s in self.products if s.price]) * self.operating_hours
    @property
    def material_cost(self):
        """Annual material cost."""
        return sum([s.cost for s in self.feeds if s.price]) * self.operating_hours
    @property
    def utility_cost(self):
        """Total utility cost in USD/yr."""
        return sum([u.utility_cost for u in self.cost_units]) * self.operating_hours
    @property
    def purchase_cost(self):
        """Total purchase cost in USD."""
        return sum([u.purchase_cost for u in self.cost_units])
    @property
    def installed_equipment_cost(self):
        """Total installed cost (USD)."""
        lang_factor = self.lang_factor
        if lang_factor:
            return sum([u.purchase_cost * lang_factor for u in self.cost_units])
        else:
            return sum([u.installed_cost for u in self.cost_units])

    def get_electricity_consumption(self):
        """Return the total electricity consumption in kWhr / yr."""
        return self.operating_hours * self.power_utility.consumption

    def get_electricity_production(self):
        """Return the total electricity production in kWhr / yr."""
        return self.operating_hours * self.power_utility.production
    
    def get_utility_duty(self, agent):
        """Return the total utility duty for given agent in kJ/yr."""
        if not isinstance(agent, str): agent = agent.ID
        return self.operating_hours * sum([i.duty for i in self.heat_utilities if i.agent.ID == agent]) 
    
    def get_utility_flow(self, agent):
        """Return the total utility flow for given agent in kJ/yr."""
        if not isinstance(agent, str): agent = agent.ID
        return self.operating_hours * sum([i.flow for i in self.heat_utilities if i.agent.ID == agent]) 
    
    def get_cooling_duty(self):
        """Return the total cooling duty in kJ/yr."""
        return - self.operating_hours * sum([i.duty for i in self.heat_utilities if i.flow * i.duty < 0])
    
    def get_heating_duty(self):
        """Return the total heating duty in kJ/yr."""
        return self.operating_hours * sum([i.duty for i in self.heat_utilities if i.flow * i.duty > 0])
    
    # Other
    def to_network(self):
        """Return network that defines the system path."""
        isa = isinstance
        path = [(i.to_network() if isa(i, System) else i) for i in self._path]
        network = Network.__new__(Network)
        network.path = path
        network.recycle = self._recycle
        network.units = set(self.unit_path)
        return network

    # Debugging
    def _turn_on(self, mode):
        """Turn on special simulation modes like `profile` or `debug`."""
        if not isinstance(mode, str):
            raise TypeError(f"mode must be a string; not a {type(mode).__name__} object")
        mode = mode.lower()
        if mode == 'debug':
            _wrap_method = _method_debug
        elif mode == 'profile':
            _wrap_method = _method_profile
        else:
            raise ValueError(f"mode must be either 'debug' or 'profile'; not '{mode}'")
        for u in self.units:
            if u._specification:
                u._specification = [_wrap_method(u, i) for i in u.specification]
            else:
                u.run = _wrap_method(u, u.run)
            u._design = _wrap_method(u, u._design)
            u._cost = _wrap_method(u, u._cost)

    def _turn_off(self):
        """Turn off special simulation modes like `profile` or `debug`."""
        for u in self.units:
            if u.specification:
                u.specification = u.specification._original
            else:
                u.run = u.run._original
            u._design = u._design._original
            u._cost = u._cost._original

    def debug(self):
        """Simulate in debug mode. If an exception is raised, it will
        automatically enter in a breakpoint"""
        self._turn_on('debug')
        try: self.simulate()
        finally: self._turn_off()

    def profile(self):
        """
        Simulate system in profile mode and return a DataFrame object of unit
        operation simulation times.

        """
        import pandas as pd
        self._turn_on('profile')
        try: self.simulate()
        finally: self._turn_off()
        units = self.units
        units.sort(key=(lambda u: u._total_excecution_time_), reverse=True)
        data = [(u.line, 1000. * u._total_excecution_time_) for u in units]
        for u in units: del u._total_excecution_time_
        return pd.DataFrame(data, index=[u.ID for u in units],
                            columns=('Unit Operation', 'Time (ms)'))

    # Representation
    def print(self, spaces=''): # pragma: no cover
        """
        Print in a format that you can use recreate the system.
        """
        print(self._stacked_info())

    def _stacked_info(self, spaces=''): # pragma: no cover
        """
        Return info with inner layers of path and facilities stacked.
        """
        info = f"{type(self).__name__}({repr(self.ID)}"
        spaces += 4 * " "
        dlim = ',\n' + spaces
        update_info = lambda new_info: dlim.join([info, new_info])
        def get_path_info(path):
            isa = isinstance
            path_info = []
            for i in path:
                if isa(i, Unit):
                    path_info.append(str(i))
                elif isa(i, System):
                    path_info.append(i._stacked_info(spaces))
                else:
                    path_info.append(str(i))
            return '[' + (dlim + " ").join(path_info) + ']'
        path_info = get_path_info(self._path)
        info = update_info(path_info)
        facilities = self._facilities
        if facilities:
            facilities_info = get_path_info(facilities)
            facilities_info = f'facilities={facilities_info}'
            info = update_info(facilities_info)
        recycle = self._recycle
        if recycle:
            recycle = self._get_recycle_info()
            info = update_info(f"recycle={recycle}")
        if self.N_runs:
            info = update_info(f"N_runs={self.N_runs}")
        info += ')'
        return info

    def _get_recycle_info(self):
        recycle = self._recycle
        if isinstance(recycle, Stream):
            recycle = recycle._source_info()
        else:
            recycle = ", ".join([i._source_info() for i in recycle])
            recycle = '{' + recycle + '}'
        return recycle

    def _ipython_display_(self):
        if bst.ALWAYS_DISPLAY_DIAGRAMS: self.diagram('minimal')
        self.show()

    def _error_info(self):
        """Return information on convergence."""
        recycle = self._recycle
        if recycle:
            s = '' if isinstance(recycle, Stream) else 's'
            return (f"\nHighest convergence error among components in recycle"
                    f"\nstream{s} {self._get_recycle_info()} after {self._iter} loops:"
                    f"\n- flow rate   {self._mol_error:.2e} kmol/hr ({self._rmol_error*100.:.2g}%)"
                    f"\n- temperature {self._T_error:.2e} K ({self._rT_error*100.:.2g}%)")
        else:
            return ""

    def __str__(self):
        if self.ID: return self.ID
        else: return type(self).__name__

    def __repr__(self):
        if self.ID: return f'<{type(self).__name__}: {self.ID}>'
        else: return f'<{type(self).__name__}>'

    def show(self, layout=None, T=None, P=None, flow=None, composition=None, N=None,
             IDs=None, data=True):
        """Prints information on system."""
        print(self._info(layout, T, P, flow, composition, N, IDs, data))

    def _info(self, layout, T, P, flow, composition, N, IDs, data):
        """Return string with all specifications."""
        error = self._error_info()
        ins_and_outs = repr_ins_and_outs(layout, self.ins, self.outs,
                                         T, P, flow, composition, N, IDs, data)
        return (f"System: {self.ID}"
                + error + '\n'
                + ins_and_outs)

class FacilityLoop(System):
    __slots__ = ()

    def _run(self):
        obj = super()
        for i in self.units:
            if i._design or i._cost: Unit._setup(i)
        obj._run()
        self._summary()

from biosteam import _flowsheet as flowsheet_module
del ignore_docking_warnings


# %% Working with different operation modes

class OperationModeResults:
  
    __slots__ = ('unit_capital_costs', 
                 'utility_cost', 
                 'stream_properties',
                 'feeds', 'products',
                 'heat_utilities',
                 'power_utility')
    
    def __init__(self, unit_capital_costs, stream_properties, 
                 utility_cost, feeds, products, heat_utilities, power_utility):
        self.unit_capital_costs = unit_capital_costs
        self.stream_properties = stream_properties
        self.utility_cost = utility_cost
        self.feeds = feeds
        self.products = products
        self.heat_utilities = heat_utilities
        self.power_utility = power_utility

    @property
    def material_cost(self):
        flow_rates = self.stream_properties['F_mass']
        return sum([flow_rates[i] * i.price  for i in self.feeds])

    @property
    def sales(self):
        flow_rates = self.stream_properties['F_mass']
        return sum([flow_rates[i] * i.price for i in self.products])


class OperationMode:
    __slots__ = ('__dict__',)
    def __init__(self, **data):
        self.__dict__ = data

    def simulate(self):
        """
        Simulate operation mode and return an OperationModeResults object with
        data on variable operating costs (i.e. utility and material costs) and sales.

        """
        agile_system = self.agile_system
        operation_parameters = agile_system.operation_parameters
        mode_operation_parameters = agile_system.mode_operation_parameters
        stream_properties = agile_system.stream_properties
        for name, value in self.__dict__.items():
            if name in operation_parameters: operation_parameters[name](value)
            elif name in mode_operation_parameters: mode_operation_parameters[name](value, self)
        system = self.system
        system.simulate()
        feeds = system.feeds
        products = system.products
        cost_units = system.cost_units
        operating_hours = self.operating_hours
        streams = feeds + products
        return OperationModeResults(
            {i: i.get_design_and_capital() for i in cost_units},
            {name: {stream: getattr(stream, name) for stream in streams}
             for name in stream_properties},
            operating_hours * sum([i.utility_cost for i in cost_units]),
            feeds, products, system.heat_utilities, system.power_utility
        )

    def __repr__(self):
        return f"{type(self).__name__}({repr_kwargs(self.__dict__, start='')})"



class OperationMetric:
    __slots__ = ('getter', 'value')
    
    def __init__(self, getter):
        self.getter = getter
    
    def __call__(self):
        return self.value
    
    def __repr__(self):
        return f"{type(self).__name__}(getter={self.getter})"



class AgileSystem:
    """
    Class for creating objects which may serve to retrive
    general results from multiple operation modes in such a way that it
    represents an agile production process. When simulated, an AgileSystem
    generates results from system operation modes and compile them to
    retrieve results later.

    Parameters
    ----------
    operation_modes : list[OperationMode], optional
        Defines each mode of operation with time steps and parameter values
    operation_parameters : dict[str: function], optional
        Defines all parameters available for all operation modes.
    lang_factor : float, optional
        Lang factor for getting fixed capital investment from
        total purchase cost. If no lang factor, installed equipment costs are
        estimated using bare module factors.

    """

    __slots__ = (
        'operation_modes', 'operation_parameters',
        'mode_operation_parameters', 'annual_operation_metrics',
        'operation_metrics', 'unit_capital_costs', 
        'net_electricity_consumption', 'utility_cost', 
        'stream_properties', 'flow_rates', 'feeds', 'products', 'purchase_cost', 
        'installed_equipment_cost', 'heat_utilities', 'power_utility',
        'process_impact_items', 'lang_factor', '_OperationMode', 
        '_TEA', '_LCA', '_units', '_streams',
    )

    TEA = System.TEA
    LCA = System.LCA
 
    define_process_impact = System.define_process_impact
    get_net_heat_utility_impact = System.get_net_heat_utility_impact
    get_net_electricity_impact = System.get_net_electricity_impact
    get_net_utility_impact = System.get_net_utility_impact
    get_total_input_impact = System.get_total_input_impact
    get_process_impact = System.get_process_impact
    get_net_impact = System.get_net_impact
    get_property_allocated_impact = System.get_property_allocated_impact
    get_property_allocation_factors = System.get_property_allocation_factors
    get_displacement_allocation_factors = System.get_displacement_allocation_factors
    get_electricity_consumption = System.get_electricity_consumption
    get_electricity_production = System.get_electricity_production
    get_utility_duty = System.get_utility_duty
    get_utility_flow = System.get_utility_flow
    get_cooling_duty = System.get_cooling_duty
    get_heating_duty = System.get_heating_duty

    def __init__(self, operation_modes=None, operation_parameters=None, 
                 mode_operation_parameters=None, annual_operation_metrics=None,
                 operation_metrics=None, lang_factor=None, 
                 stream_property_names=None):
        self.operation_modes = [] if operation_modes is None else operation_modes 
        self.operation_parameters = {} if operation_parameters  is None else operation_parameters
        self.mode_operation_parameters = {} if mode_operation_parameters is None else mode_operation_parameters
        self.annual_operation_metrics = [] if annual_operation_metrics is None else annual_operation_metrics
        self.operation_metrics = [] if operation_metrics is None else operation_metrics
        self.flow_rates = flow_rates = {}
        self.stream_properties = stream_properties = {'F_mass': flow_rates}
        if stream_property_names is not None:
            for i in stream_property_names: stream_properties[i] = {}
        self.lang_factor = lang_factor
        self.heat_utilities = None
        self.power_utility = None
        self._OperationMode = type('OperationMode', (OperationMode,), {'agile_system': self})

    def _downstream_system(self, unit):
        return self

    def operation_mode(self, system, operating_hours, **data):
        """
        Define and register an operation mode.

        Parameters
        ----------
        operating_hours : function
            Length of operation in hours.
        **data : str
            Name and value-pairs of operation parameters.

        """
        for s in system.streams: s.unlink()
        om = self._OperationMode(system=system, operating_hours=operating_hours, **data)
        self.operation_modes.append(om)
        return om

    def operation_parameter(self, setter=None, name=None, mode_dependent=False):
        """
        Define and register operation parameter.

        Parameters
        ----------
        setter : function
            Should set parameter.
        name : str, optional
            Name of parameter. If None, default to argument name of setter.
        mode_dependent :
            Whether the setter accepts the OperationMode object as a second argument.

        """
        if not setter: return lambda setter: self.operation_parameter(setter, name, mode_dependent)
        if not name: name, *_ = signature(setter).parameters.keys()
        if mode_dependent:
            self.mode_operation_parameters[name] = setter
        else:
            self.operation_parameters[name] = setter
        return setter

    def operation_metric(self, getter=None, annualize=False):
        """
        Return an OperationMetric object.
        
        Parameters
        ----------    
        getter : Callable(<OperationMode>)
            Should return the metric value.
        annualize : bool, optional
            Whether to multiply by operating hours and sum across operation modes.
            If True, return value of the OperationMetric is a float. Otherwise, 
            the OperationMetric returns a dictionary of results with OperationMode
            objects as keys.
        
        """
        if not getter: return lambda getter: self.operation_metric(getter, annualize)
        operation_metric = OperationMetric(getter)
        if annualize:
            self.annual_operation_metrics.append(operation_metric)
        else:
            self.operation_metrics.append(operation_metric)
        return operation_metric

    def get_market_value(self, stream):
        """Return the market value of a stream [USD/yr]."""
        return self.flow_rates[stream] * stream.price
  
    def get_mass_flow(self, stream):
        """Return the mass flow rate of a stream [kg/yr]."""
        return self.flow_rates[stream]
    
    def get_property(self, stream, name, units=None):
        """Return the annualized property of a stream."""
        value = self.stream_properties[name][stream] * self.operating_hours
        if units is None:
            return value
        else:
            units_dct = bst.Stream._units_of_measure
            original_units = units_dct[name]
            return original_units.convert(value, units)
    
    def get_total_feeds_impact(self, key):
        """
        Return the total annual impact of all feeds given 
        the characterization factor key.
        
        """
        flow_rates = self.flow_rates
        return sum([flow_rates[s] * s.characterization_factors[key] for s in self.feeds
                    if key in s.characterization_factors])
    
    def get_total_products_impact(self, key):
        """
        Return the total annual impact of all products given 
        the characterization factor key.
        
        """
        flow_rates = self.flow_rates
        return sum([flow_rates[s] * s.characterization_factors[key] for s in self.products
                    if key in s.characterization_factors])
    
    def get_material_impact(self, stream, key):
        """
        Return the annual material impact given the stream and the 
        characterization factor key.
        
        """
        return self.flow_rates[stream] * stream.characterization_factor[key]
    
    def _price2cost(self, stream):
        """Get factor to convert stream price to cost for cash flow in solve_price method."""
        if stream in self.flow_rates:
            F_mass = self.flow_rates[stream]
        else:
            F_mass = 0.
        if not F_mass: warn(f"stream '{stream}' is empty", category=RuntimeWarning)
        if stream in self.products:
            return F_mass
        elif stream in self.feeds:
            return - F_mass
        else:
            raise ValueError("stream must be either a feed or a product")

    @property
    def material_cost(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price  for i in self.feeds])

    @property
    def sales(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price for i in self.products])

    @property
    def streams(self):
        try:
            return self._streams
        except:
            self._streams = streams = []
            stream_set = set()
            for u in self.units:
                for s in u._ins + u._outs:
                    if not s or s in stream_set: continue
                    streams.append(s)
                    stream_set.add(s)
            return streams

    @property
    def units(self):
        try:
            return self._units
        except:
            self._units = units = []
            past_units = set()
            for i in self.operation_modes:
                for i in i.system.unit_path:
                    if i in past_units: continue
                    units.append(i)
            return units

    @property
    def cost_units(self):
        systems = set([i.system for i in self.operation_modes])
        if len(systems) == 1:
            return systems.pop().cost_units
        else:
            units = set()
            for i in systems: units.update(i.cost_units)
            return units

    def empty_recycles(self):
        for mode in self.operation_modes:
            mode.system.empty_recycles()

    def reset_cache(self):
        for mode in self.operation_modes:
            mode.system.reset_cache()

    @property
    def operating_hours(self):
        return sum([i.operating_hours for i in self.operation_modes])
    @operating_hours.setter
    def operating_hours(self, operating_hours):
        factor = operating_hours / self.operating_hours
        for i in self.operation_modes: i.operating_hours *= factor

    def reduce_chemicals(self, required_chemicals=()):
        for i in self.streams: i.unlink()
        unit_thermo = {}
        mixer_thermo = {}
        thermo_cache = {}
        for mode in self.operation_modes:
            mode.system._load_configuration()
            mode.system._reduced_thermo_data(required_chemicals, unit_thermo, mixer_thermo, thermo_cache)
        for mode in self.operation_modes:
            mode.system._load_configuration()
            for unit, thermo in unit_thermo.items(): unit._reset_thermo(thermo)
            for mixer, thermo in mixer_thermo.items():
                for i in mixer._ins:
                    if i._source: i._reset_thermo(unit_thermo[i._source])
                thermo = mixer_thermo[mixer]
                mixer._load_thermo(thermo)
                mixer._outs[0]._reset_thermo(thermo)
            mode.system._interface_property_packages()

    def simulate(self):
        operation_modes = self.operation_modes
        operation_metrics = self.operation_metrics
        annual_operation_metrics = self.annual_operation_metrics
        N_modes = len(operation_modes)
        N_metrics = len(operation_metrics)
        N_annual_metrics = len(annual_operation_metrics)
        annual_metric_range = range(N_annual_metrics)
        metric_range = range(N_metrics)
        mode_range = range(N_modes)
        operation_mode_results = N_modes * [None]
        annual_values = [N_modes * [None] for i in annual_metric_range]
        values = [{i: None for i in operation_modes} for i in metric_range]
        total_operating_hours = self.operating_hours
        for i in mode_range:
            mode = operation_modes[i]
            operation_mode_results[i] = results = mode.simulate()
            for j in annual_metric_range:
                metric = annual_operation_metrics[j]
                annual_values[j][i] = metric.getter(mode) * mode.operating_hours
            for j in metric_range:
                metric = operation_metrics[j]
                values[j][mode] = metric.getter(mode)
            scale = mode.operating_hours / total_operating_hours
            for hu in results.heat_utilities: hu.scale(scale)
            results.power_utility.scale(scale)
        for i in annual_metric_range:
            metric = annual_operation_metrics[i]
            metric.value = sum(annual_values[i])
        for i in metric_range:
            metric = operation_metrics[i]
            metric.value = values[i]
        units = set(sum([list(i.unit_capital_costs) for i in operation_mode_results], []))
        unit_modes = {i: [] for i in units}
        for results in operation_mode_results:
            for i, j in results.unit_capital_costs.items(): unit_modes[i].append(j)
        self.heat_utilities = bst.HeatUtility.sum_by_agent(sum([r.heat_utilities for r in operation_mode_results], ()))
        self.power_utility = bst.PowerUtility.sum([r.power_utility for r in operation_mode_results])
        self.unit_capital_costs = {i: i.get_agile_design_and_capital(j) for i, j in unit_modes.items()}
        self.utility_cost = sum([i.utility_cost for i in operation_mode_results])
        self.feeds = list(set(sum([i.feeds for i in operation_mode_results], [])))
        self.products = list(set(sum([i.products for i in operation_mode_results], [])))
        self.purchase_cost = sum([u.purchase_cost for u in self.unit_capital_costs])
        lang_factor = self.lang_factor
        if lang_factor:
            self.installed_equipment_cost = sum([u.purchase_cost * lang_factor for u in self.unit_capital_costs.values()])
        else:
            self.installed_equipment_cost = sum([u.installed_cost for u in self.unit_capital_costs.values()])
        stream_properties = self.stream_properties
        for dct in stream_properties.values(): dct.clear()
        for i in mode_range:
            results = operation_mode_results[i]
            operating_hours = operation_modes[i].operating_hours
            for name, dct in results.stream_properties.items():
                propdct = stream_properties[name]
                for stream, property in dct.items():
                    if stream in propdct: propdct[stream] += property * operating_hours
                    else: propdct[stream] = property * operating_hours

    def __repr__(self):
        return f"{type(self).__name__}(operation_modes={self.operation_modes}, operation_parameters={self.operation_parameters}, lang_factor={self.lang_factor})"