# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
from collections.abc import Iterable
from warnings import warn
import biosteam as bst
import numpy as np

__all__ = ('System', 'MockSystem', 'ScenarioCosts', 
           'mark_disjunction', 'unmark_disjunction')    

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


# %% Functions for recycle

def check_recycle_feasibility(material: np.ndarray):
    if fn.infeasible(material):
        raise InfeasibleRegion('recycle material flow rate')
    else:
        material[material < 0.] = 0. 


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
        if hasattr(self, '_ins'):
            ins = self._ins
        else:
            inlets = bst.utils.feeds_from_units(self.units)
            self._ins = ins = StreamPorts.from_inlets(inlets, sort=True)
        return ins
    @property
    def outs(self):
        """StreamPorts[:class:`~OutletPort`] All outlets to the system."""
        if hasattr(self, '_outs'):
            outs = self._outs
        else:
            outlets = bst.utils.products_from_units(self.units)
            self._outs = outs = StreamPorts.from_outlets(outlets, sort=True)
        return outs
    
    def load_inlet_ports(self, inlets):
        """Load inlet ports to system."""
        all_inlets = bst.utils.inlets(self.units)
        for i in inlets: 
            if i not in all_inlets:
                raise ValueError(f'{i} is not an inlet')
        self._ins = StreamPorts.from_inlets(inlets)
    
    def load_outlet_ports(self, outlets):
        """Load outlet ports to system."""
        all_outlets = bst.utils.outlets(self.units)
        for i in outlets: 
            if i not in all_outlets:
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
        '_connections',
        '_irrelevant_units',
        '_converge_method',
        '_TEA',
        '_LCA',
    )
    
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

    #: [bool] Whether to raise a RuntimeError when system doesn't converge
    strict_convergence = True

    @classmethod
    def from_feedstock(cls, ID, feedstock, feeds=None, facilities=(), 
                       ends=None, facility_recycle=None, operating_hours=None):
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
        
        """
        network = Network.from_feedstock(feedstock, feeds, ends)
        return cls.from_network(ID, network, facilities, 
                                facility_recycle, operating_hours)

    @classmethod
    def from_units(cls, ID="", units=None, feeds=None, ends=None,
                   facility_recycle=None, operating_hours=None):
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
                operating_hours=operating_hours,
            )
        else:
            system = cls(ID, (), operating_hours=operating_hours)
        return system

    @classmethod
    def from_network(cls, ID, network, facilities=(), facility_recycle=None,
                     operating_hours=None):
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
        self._load_stream_links()
        self._load_defaults()
        self._save_configuration()
        if operating_hours is not None: self.operating_hours = operating_hours
        return self
        
    def __init__(self, ID, path=(), recycle=None, facilities=(), 
                 facility_recycle=None, N_runs=None, operating_hours=None):
        self.recycle = recycle
        self.N_runs = N_runs
        self._set_path(path)
        self._specification = None
        self._load_flowsheet()
        self._reset_errors()
        self._set_facilities(facilities)
        self._set_facility_recycle(facility_recycle)
        self._load_stream_links()
        self._register(ID)
        self._load_defaults()
        self._save_configuration()
        if operating_hours is not None: self.operating_hours = operating_hours
    
    def __enter__(self):
        if self.path or self.recycle or self.facilities:
            raise RuntimeError("only empty systems can enter `with` statement")
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
        if self.path or self.recycle or self.facilities:
            raise RuntimeError('system cannot be modified before exiting `with` statement')
        else:
            units = [i for i in dump if i not in irrelevant_units]
            system = self.from_units(None, units)
            self.ID = ID
            self.copy_like(system)
    
    def _save_configuration(self):
        self._connections = [i.get_connection() for i in self.streams]
    
    @ignore_docking_warnings
    def _load_configuration(self):
        for i in self._connections:
            if i.source:
                i.source.outs[i.source_index] = i.stream
            if i.sink:
                i.sink.ins[i.sink_index] = i.stream    
    
    def copy_like(self, other):
        """Copy path, facilities and recycle from other system.""" 
        self._path = other._path
        self._facilities = other._facilities
        self._facility_loop = other._facility_loop
        self._recycle = other._recycle
        self._connections = other._connections
    
    def set_tolerance(self, mol=None, rmol=None, T=None, rT=None, subsystems=False):
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

        """
        if mol: self.molar_tolerance = float(mol)
        if rmol: self.relative_molar_tolerance = float(rmol)
        if T: self.temperature_tolerance = float(T)
        if rT: self.temperature_tolerance = float(rT)
        if subsystems: 
            for i in self.subsystems: i.set_tolerance(mol, rmol, T, rT, subsystems)
    
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
    
    specification = Unit.specification
    save_report = save_report
    
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
        if self.recycle: raise RuntimeError('cannot split system with recycle')
        path = self.path
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
                System(ID_downstream, path[index:], None, self.facilities))
    
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
                i._system = self
                i._other_units = other_units = units.copy()
                other_units.remove(i)
            
            
    def _set_facility_recycle(self, recycle):
        if recycle:
            sys = self._downstream_system(recycle.sink)
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
        return [i for i in self.path if isinstance(i, System)]
    @property
    def units(self):
        """[list] All unit operations as ordered in the path."""
        units = []
        isa = isinstance
        for i in self._path + self._facilities:
            if isa(i, Unit):
                units.append(i)
            elif isa(i, System):
                units.extend(i.units)
        return units 
    @property
    def streams(self):
        """set[:class:`~thermosteam.Stream`] All streams within the system."""
        streams = bst.utils.streams_from_units(self.units)
        bst.utils.filter_out_missing_streams(streams)
        return streams
    @property
    def feeds(self):
        """set[:class:`~thermosteam.Stream`] All feeds to the system."""
        inlets = bst.utils.inlets(self.units)
        bst.utils.filter_out_missing_streams(inlets)
        return bst.utils.feeds(inlets)
    @property
    def products(self):
        """set[:class:`~thermosteam.Stream`] All products of the system."""
        outlets = bst.utils.outlets(self.units)
        bst.utils.filter_out_missing_streams(outlets)
        return bst.utils.products(outlets)
    
    @property
    def TEA(self):
        """[TEA] Object for Techno-Economic Analysis."""
        try: return self._TEA
        except AttributeError: return None
    
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

    def _downstream_path(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.units: return []
        elif self._recycle: return self._path
        isa = isinstance
        for index, i in enumerate(self._path):
            if unit is i:
                return self._path[index:]
            elif (isa(i, System) and unit in i.units): 
                return i._downstream_path(unit) + self._path[index+1:]
        return []
    
    def _downstream_facilities(self, unit):
        """Return a list of facilities composed of the `unit` and 
        everything downstream."""
        isa = isinstance
        for index, i in enumerate(self._facilities):
            if unit is i or (isa(i, System) and unit in i.units):
                return self._facilities[index:]
        return []
    
    def _downstream_system(self, unit):
        """Return a system with a path composed of the `unit` and
        everything downstream (facilities included)."""
        if self.recycle or unit is self._path[0]: return self
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
        return digraph_from_units_and_streams(self.units, 
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
        check_recycle_feasibility(mol)
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
                self._rmol_error = rmol_error = (mol_errors / np.maximum.reduce([mol[positive_index], mol_new[positive_index]])).max()
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
            return recycle.imol.data.copy()
        elif isinstance(recycle, Iterable):
            return np.vstack([i.imol.data for i in recycle])
        else:
            raise RuntimeError('no recycle available')
    
    def _set_recycle_data(self, data):
        recycle = self._recycle
        isa = isinstance
        if isa(recycle, Stream):
            try:
                recycle._imol._data[:] = data
            except:
                raise IndexError(f'expected 1 row; got {data.shape[0]} rows instead')
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
        try:
            solver(self._iter_run, self._get_recycle_data())
        except IndexError as error:
            try:
                solver(self._iter_run, self._get_recycle_data())
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

    def reset_cache(self):
        """Reset cache of all unit operations."""
        for unit in self.units: unit.reset_cache()

    def simulate(self):
        """Converge the path and simulate all units."""
        self._setup()
        self._converge()
        self._summary()
        if self._facility_loop: self._facility_loop._converge()
    
    # Convinience methods
    
    @property
    def heat_utilities(self):
        """[tuple] All HeatUtility objects."""
        return utils.get_heat_utilities(self.units)
    
    @property
    def power_utilities(self):
        """[tuple] All PowerUtility objects."""
        return tuple(utils.get_power_utilities(self.units))
    
    def get_inlet_flow(self, units, key=None):
        """
        Return total flow across all inlets.
        
        Parameters
        ----------
        units : str
            Units of measure.
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned
            
        Examples
        --------
        >>> from biorefineries.cornstover import cornstover_sys
        >>> from biosteam import default
        >>> cornstover_sys.get_inlet_flow('tonne/s') # Sum of all chemicals
        51422.13
        >>> cornstover_sys.get_inlet_flow('tonne/s', 'Water') # Just water
        46050.96
        >>> default() # Bring biosteam settings back to default
        
        """
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in bst.utils.inlets(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in bst.utils.inlets(self.units)])
    
    def get_outlet_flow(self, units, key=None):
        """
        Return total flow across all outlets.
        
        Parameters
        ----------
        units : str
            Units of measure.
        key : tuple[str] or str, optional
            Chemical identifiers. If none given, the sum of all chemicals returned
            
        Examples
        --------
        >>> from biorefineries.cornstover import cornstover_sys
        >>> from biosteam import default
        >>> cornstover_sys.get_outlet_flow('tonne/s') # Sum of all chemicals
        51558.88
        >>> cornstover_sys.get_outlet_flow('tonne/s', 'Water') # Just water
        46103.38
        >>> default() # Bring biosteam settings back to default
        
        """
        if key:
            return self.operating_hours * sum([i.get_flow(units, key) for i in bst.utils.outlets(self.units)])
        else:
            return self.operating_hours * sum([i.get_total_flow(units) for i in bst.utils.outlets(self.units)])
    
    def get_electricity_consumption(self):
        """Return the total electricity consumption in MW."""
        return self.operating_hours * utils.get_electricity_consumption(self.power_utilities)

    def get_electricity_production(self):
        """Return the total electricity production in MW."""
        return self.operating_hours * utils.get_electricity_production(self.power_utilities)
    
    def get_utility_duty(self, agent):
        """Return the total utility duty for given agent in GJ/hr"""
        return self.operating_hours * utils.get_utility_duty(self.heat_utilities, agent)
    
    def get_utility_flow(self, agent):
        """Return the total utility flow for given agent in MT/hr"""
        return self.operating_hours * utils.get_utility_flow(self.heat_utilities, agent)
    
    def get_cooling_duty(self):
        """Return the total cooling duty in GJ/hr."""
        return self.operating_hours * utils.get_cooling_duty(self.heat_utilities)
    
    def get_heating_duty(self):
        """Return the total heating duty in GJ/hr."""
        return self.operating_hours * utils.get_heating_duty(self.heat_utilities)
    
    def get_purchase_cost(self):
        """Return the total equipment purchase cost in million USD."""
        return utils.get_purchase_cost(self.units)
    
    def get_installed_equipment_cost(self):
        """Return the total installed equipment cost in million USD."""
        return utils.get_installed_cost(self.units)
    
    def get_scenario_costs(self):
        """
        Return a ScenarioCosts object with data on variable operating costs
        (i.e. utility and material costs) and sales.
        
        Examples
        --------
        Create a simple heat exchanger system and get scenario costs:
            
        >>> import biosteam as bst
        >>> bst.default() # Reset to biosteam defaults
        >>> bst.settings.set_thermo(['Water'])
        >>> feed = bst.Stream('feed', Water=100)
        >>> product = bst.Stream('product')
        >>> HX1 = bst.HXutility('HX1', feed, product, T=350)
        >>> operating_days = 300
        >>> operating_hours = 24 * operating_days
        >>> SYS1 = bst.System.from_units('SYS1', [HX1], operating_hours=operating_hours)
        >>> SYS1.simulate()
        >>> scenario_costs = SYS1.get_scenario_costs()
        
        A ScenarioCosts object has useful properties that may change according to
        stream prices, but not with flow rates or unit simulations:
            
        >>> assert scenario_costs.utility_cost == HX1.utility_cost * operating_hours
        >>> assert scenario_costs.material_cost == feed.cost == 0.
        >>> assert scenario_costs.sales == product.cost == 0.
        >>> # Note how the feed price affects material costs
        >>> feed.price = 0.05
        >>> assert scenario_costs.material_cost == feed.cost * operating_hours
        >>> # Yet simulations and changes to material flow rates do not affect costs
        >>> feed.set_total_flow(200, 'kmol/hr')
        >>> HX1.T = 370
        >>> HX1.simulate()
        >>> assert scenario_costs.material_cost != feed.cost * operating_hours
        >>> assert scenario_costs.utility_cost != HX1.utility_cost * operating_hours
        
        """
        feeds = self.feeds
        products = self.products
        units = self.units
        operating_hours = self.operating_hours
        return ScenarioCosts(
            {i: i.get_capital_costs() for i in units if i._design or i._cost},
            {i: i.F_mass * operating_hours for i in feeds + products},
            operating_hours * sum([i.utility_cost for i in units]),
            feeds, products,
            self.operating_hours,
        )
    
    # Other
    def to_network(self):
        """Return network that defines the system path."""
        isa = isinstance
        path = [(i.to_network() if isa(i, System) else i) for i in self._path]
        network = Network.__new__(Network)    
        network.path = path
        network.recycle = self._recycle
        network.units = set(self.units)
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
            if u.specification:
                u.specification = _wrap_method(u, u.specification)
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
        facilities = self.facilities
        if facilities:
            facilities_info = get_path_info(self.facilities)
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


# %% Working with different scenarios

class ScenarioCosts:
    
    __slots__ = ('unit_capital_costs', 'utility_cost', 'flow_rates', 
                 'feeds', 'products', 'operating_hours')
    
    def __init__(self, unit_capital_costs, flow_rates, utility_cost, 
                 feeds, products, operating_hours):
        self.unit_capital_costs = unit_capital_costs
        self.flow_rates = flow_rates
        self.utility_cost = utility_cost
        self.feeds = feeds
        self.products = products
        self.operating_hours = operating_hours
    
    @property
    def material_cost(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price  for i in self.feeds])
    
    @property
    def sales(self):
        flow_rates = self.flow_rates
        return sum([flow_rates[i] * i.price for i in self.products])