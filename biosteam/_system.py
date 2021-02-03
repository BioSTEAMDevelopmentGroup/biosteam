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
from .report import save_report
from .exceptions import InfeasibleRegion
from .utils import StreamPorts, OutletPort, colors
from collections.abc import Iterable
from collections import deque
from warnings import warn
import biosteam as bst
import numpy as np

__all__ = ('System', 'MockSystem', 'mark_disjunction', 'unmark_disjunction')    

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
    
def run_unit_in_path(unit):
    specification = unit._specification
    if specification:
        method = specification
    else:
        method = unit._run
    try_method_with_object_stamp(unit, method)

def converge_system_in_path(system):
    specification = system._specification
    if specification:
        method = specification
    else:
        method = system._converge
    try_method_with_object_stamp(system, method)

def simulate_unit_in_path(unit):
    try_method_with_object_stamp(unit, unit.simulate)

def simulate_system_in_path(system):
    specification = system._specification
    if specification:
        method = specification
    else:
        method = system.simulate
    try_method_with_object_stamp(system, method)

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
    f = bst.main_flowsheet
    for i in (f.system, f.stream, f.unit, f.flowsheet): lcs.update(i.__dict__)
    lcs.update(bst.__dict__)

def _method_debug(self, func):
    """Method decorator for debugging system."""
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            print_exception_in_debugger(self, func, e)
            update_locals_with_flowsheet(locals())
            
            # Greatings from the Yoel, the BDFL of BioSTEAM.
            # All systems, units, streams, and flowsheets are available as 
            # local variables. Although this debugging method is meant to 
            # be for internal development, please feel free to give it a shot.
            breakpoint()
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


# %% Converging recycle systems

class MockSystem:
    """
    Create a MockSystem object with inlets and outlets just like System 
    objects, but without implementing any of the convergence methods nor
    path related attributes.
    
    Notes
    -----
    This object is used to prevent the creation of unneeded systems for less 
    computational effort.
    
    """
    __slots__ = ('_ins', '_outs')
    
    def __init__(self, ins, outs):
        self._ins = StreamPorts.from_inlets(ins)
        self._outs = StreamPorts.from_outlets(outs)
        
    @property
    def ins(self): return self._ins
    @property
    def outs(self): return self._outs
        
    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    def show(self):
        print(
            f"{type(self).__name__}(\n"
            f"    ins={self.ins},\n"
            f"    outs={self.outs}\n"
             ")"
        )
        
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
    path : tuple[Unit, function and/or System], optional
        A path that is run element by element until the recycle converges.
    recycle=None : :class:`~thermosteam.Stream` or tuple[:class:`~thermosteam.Stream`], optional
        A tear stream for the recycle loop.
    facilities=() : tuple[Unit, function, and/or System], optional
        Offsite facilities that are simulated only after
        completing the path simulation.
    facility_recycle : :class:`~thermosteam.Stream`, optional
        Recycle stream between facilities and system path.
    N_runs : int, optional
        Number of iterations to run system. This parameter is applicable 
        only to systems with no recycle loop.

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
        'units',
        'subsystems',
        'maxiter',
        'molar_tolerance',
        'relative_molar_tolerance',
        'temperature_tolerance',
        'relative_temperature_tolerance',
        'flowsheet',
        'alternative_convergence_check',
        '__previous_flowsheet_units',
        '_converge_method',
        '_costunits',
        '_TEA',
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

    #: [bool] Whether to allow systems to rotate for more robust convergence.
    allow_system_rotation = False

    @classmethod
    def from_feedstock(cls, ID, feedstock, feeds=None, facilities=(), 
                       ends=None, facility_recycle=None):
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
        
        """
        network = Network.from_feedstock(feedstock, feeds, 
                                         ends or [i.get_stream() for i in disjunctions])
        return cls.from_network(ID, network, facilities, facility_recycle)

    @classmethod
    def from_units(cls, ID="", units=None, feeds=None, ends=None, facility_recycle=None):
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
        
        """
        if units is None: 
            units = ()
        elif feeds is None:
            feeds = bst.utils.feeds_from_units(units)
            bst.utils.sort_feeds_big_to_small(feeds)
        if feeds:
            feedstock, *feeds = feeds
            facilities = facilities_from_units(units) if units else ()
            system = cls.from_feedstock(
                ID, feedstock, feeds, facilities, ends,
                facility_recycle or find_blowdown_recycle(facilities)
            )
        else:
            system = cls(ID, ())
        return system

    @classmethod
    def from_network(cls, ID, network, facilities=(), facility_recycle=None):
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
        
        """
        facilities = Facility.ordered_facilities(facilities)
        isa = isinstance
        name = ID if ID is None else ''
        path = [(cls.from_network(name, i) if isa(i, Network) else i)
                for i in network.path]
        self = cls.__new__(cls)
        self.units = network.units
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
        return self
        
    def __init__(self, ID, path=(), recycle=None, facilities=(), 
                 facility_recycle=None, N_runs=None):
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
    
    def __enter__(self):
        if self.path or self.recycle or self.facilities:
            raise RuntimeError("only empty systems can enter `with` statement")
        self.__previous_flowsheet_units = set(self.flowsheet.unit)
        return self
    
    def __exit__(self, type, value, traceback):
        if value: raise value
        previous_flowsheet_units = self.__previous_flowsheet_units
        ID = self._ID
        del self.__previous_flowsheet_units
        if self.path or self.recycle or self.facilities:
            system = self.flowsheet.system[ID]
            if (system is not self 
                and system.path == self.path
                and system.recycle == self.recycle
                and system.facilities == self.facilities):
                system._ID = None
                self.ID = ID
            else:
                raise RuntimeError('system was modified before exiting `with` statement')
        else:
            units = [i for i in self.flowsheet.unit
                     if i not in previous_flowsheet_units]
            system = self.from_units('', units)
            self.ID = ID
            self.copy_like(system)
    
    def copy_like(self, other):
        """Copy path, facilities and recycle from other system.""" 
        self._path = other._path
        self._facilities = other._facilities
        self._facility_loop = other._facility_loop
        self._recycle = other._recycle
        self.units = other.units
        self.subsystems = other.subsystems
        self._costunits = other._costunits
    
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
    
    def __eq__(self, other):
        return (isinstance(other, System) 
                and self._path == other._path 
                and self._facilities == other._facilities 
                and self._recycle == other._recycle)
    
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
        
        #: [function(recycle, mol, mol_new, T, T_new) -> bool] Function that returns 
        #: whether the system has not converged (and needs to keep running).
        self.alternative_convergence_check = None
        
    
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
    
    def _get_recycle_streams(self):
        recycles = []
        recycle = self._recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycles.append(recycle)
            elif isinstance(recycle, Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_error(recycle)
        for i in self.subsystems:
            recycles.extend(i._get_recycle_streams())
        return recycles
    
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
    
    def _load_flowsheet(self):
        self.flowsheet = flowsheet_module.main_flowsheet.get_flowsheet()
    
    def _set_path(self, path):
        #: tuple[Unit, function and/or System] A path that is run element
        #: by element until the recycle converges.
        self._path = path = tuple(path)
        
        #: set[Unit] All units within the system
        self.units = units = set()
        
        #: list[System] All subsystems in the system
        self.subsystems = subsystems = []
        
        #: set[Unit] All units that have costs (including facilities).
        self._costunits = costunits = set()
        
        isa = isinstance
        for i in path:
            if isa(i, Unit): 
                units.add(i)
                if i._design or i._cost:
                    costunits.add(i)
            elif isa(i, System):
                subsystems.append(i)
                units.update(i.units)
                costunits.update(i._costunits)
    
    def _set_facilities(self, facilities):
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the path simulation.
        self._facilities = facilities = tuple(facilities)
        subsystems = self.subsystems
        costunits = self._costunits
        units = self.units
        isa = isinstance
        new_facility_units = []
        for i in facilities:
            if isa(i, Unit):
                units.add(i)
                if i._cost: costunits.add(i)
                if isa(i, Facility) and not i._system:
                    new_facility_units.append(i)
            elif isa(i, System):
                units.update(i.units)
                subsystems.append(i)
                costunits.update(i._costunits)
        for i in new_facility_units:
            i._system = self
            i._other_units = other_units = self.units.copy()
            other_units.remove(i)
            
    def _set_facility_recycle(self, recycle):
        if recycle:
            system = self._downstream_system(recycle.sink)
            #: [FacilityLoop] Recycle loop for converging facilities
            self._facility_loop = FacilityLoop(system, recycle)
        else:
            self._facility_loop = None
    
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
            
    def _load_stream_links(self):
        for u in self.unit_path: u._load_stream_links()
        
    # Forward pipping
    __sub__ = Unit.__sub__
    __rsub__ = Unit.__rsub__

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
        
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
        if unit is self._path[0]: return self
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
    
    def _minimal_digraph(self, **graph_attrs):
        """Return digraph of the path as a box."""
        return minimal_digraph(self.ID, self.units, self.streams, **graph_attrs)

    def _surface_digraph(self, **graph_attrs):
        return surface_digraph(self._path)

    def _thorough_digraph(self, **graph_attrs):
        return digraph_from_units_and_streams(self.unit_path, 
                                              self.streams, 
                                              **graph_attrs)
        
    def _cluster_digraph(self, **graph_attrs):
        return digraph_from_system(self, **graph_attrs)
        
    def diagram(self, kind=None, file=None, format=None, display=True,
                **graph_attrs):
        """
        Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of 
        the system.
        
        Parameters
        ----------
        kind : 'cluster', 'thorough', 'surface', or 'minimal'
            * **'cluster':** Display all units clustered by system.
            * **'thorough':** Display every unit within the path.
            * **'surface':** Display only elements listed in the path.
            * **'minimal':** Display path as a box.
        file=None : str, display in console by default
            File name to save diagram.
        format='png' : str
            File format (e.g. "png", "svg").
        display : bool, optional
            Whether to display diagram in console or to return the graphviz 
            object.
            
        """
        if not kind: kind = 'surface'
        if not format: format = 'png'
        if kind == 'cluster':
            f = self._cluster_digraph(format=format, **graph_attrs)
        elif kind == 'thorough':
            f = self._thorough_digraph(format=format, **graph_attrs)
        elif kind == 'surface':
            f = self._surface_digraph(format=format, **graph_attrs)
        elif kind == 'minimal':
            f = self._minimal_digraph(format=format, **graph_attrs)
        else:
            raise ValueError("kind must be either 'cluster', 'thorough', 'surface', or 'minimal'")
        if display or file: 
            finalize_digraph(f, file, format)
        else:
            return f
            
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
        molar_tolerance = self.molar_tolerance
        mol_errors = np.abs(mol - mol_new)
        positive_index = mol_errors > 1e-16
        mol_errors = mol_errors[positive_index]
        if mol_errors.size == 0:
            self._mol_error = mol_error = 0.
            self._rmol_error = rmol_error = 0.
        else:
            self._mol_error = mol_error = mol_errors.max()
            self._rmol_error = rmol_error = (mol_errors / np.maximum.reduce([mol[positive_index], mol_new[positive_index]])).max()
        T_errors = np.abs(T - T_new)
        self._T_error = T_error = T_errors.max()
        self._rT_error = rT_error = (T_errors / T).max()
        self._iter += 1
        if self.alternative_convergence_check:
            not_converged = self.alternative_convergence_check(self.recycle, mol, mol_new, T, T_new)
        elif (mol_error < molar_tolerance
            and rmol_error < self.relative_molar_tolerance
            and T_error < self.temperature_tolerance
            and rT_error < self.relative_temperature_tolerance):
            not_converged = False
        else:
            not_converged = True
        if not_converged and self._iter >= self.maxiter:
            raise RuntimeError(f'{repr(self)} could not converge' + self._error_info())
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
            recycle.imol.data[:] = data
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
                    i.imol.data[:] = data[index:next_index, :]
                    index = next_index
                else:
                    i.imol.data[:] = data[index, :]
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
        isa = isinstance
        types = (Unit, System)
        for i in self._path:
            if isa(i, types): i._setup()
        
    def _run(self):
        """Run each element in the path. Rotate path if necessary."""
        if self.allow_system_rotation and self._recycle:
            N_elements = len(self._path)
            if N_elements <= 1:
                self._run_path()
            else:
                max_rotations = N_elements - 1
                for n in range(N_elements):
                    try:
                        self._run_path()
                    except Exception as error:
                        if n == max_rotations: raise error
                        if not n: self._path = deque(self._path)
                        self._path.rotate()
                    else:
                        break
                if n: self._path = tuple(self._path)
        else:
            self._run_path()
        
    def _run_path(self):
        """Rigorously run each element in the path."""
        isa = isinstance
        run = run_unit_in_path
        converge = converge_system_in_path
        for i in self._path:
            if isa(i, Unit): run(i)
            elif isa(i, System): converge(i)
            else: i() # Assume it is a function
    
    @property
    def unit_path(self):
        """[list] All unit operations as ordered in the path."""
        unit_path = []
        isa = isinstance
        for i in self._path + self._facilities:
            if isa(i, Unit):
                unit_path.append(i)
            elif isa(i, System):
                unit_path.extend(i.unit_path)
        return unit_path 
    
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
        iscallable = callable
        for i in self._path:
            if not iscallable(i): try_method_with_object_stamp(i, i._summary)
        isa = isinstance
        simulate_unit = simulate_unit_in_path
        simulate_system = simulate_system_in_path
        for i in self._facilities:
            if isa(i, Unit): simulate_unit(i)
            elif isa(i, System): simulate_system(i)
            else: i() # Assume it is a function

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
    
    def reset_flows(self):
        """Reset all process streams to zero flow."""
        from warnings import warn
        warn(DeprecationWarning("'reset_flows' will be depracated; please use 'empty_process_streams'"))
        self.empty_process_streams()
        
    def empty_process_streams(self):
        """Reset all process streams to zero flow."""
        self._reset_errors()        
        feeds = self.feeds
        for stream in self.streams:
            if stream not in feeds: stream.empty()

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
        if self._facility_loop: self._facility_loop()
     
    # Other
    
    def to_network(self):
        """Return network that defines the system path."""
        isa = isinstance
        path = [(i.to_network() if isa(i, System) else i) for i in self._path]
        network = Network.__new__(Network)    
        network.path = path
        network.recycle = self._recycle
        network.units = self.units
        network.subnetworks = [i for i in path if isa(i, Network)]
        network.feeds = self.feeds
        network.products = self.products
        return network
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        path = self._path
        for i, item in enumerate(path):
            if isinstance(item, Unit):
                item._run = _method_debug(item, item._run)
            elif callable(item):
                path[i] = _method_debug(item, item)

    def _debug_off(self):
        """Turn off debug mode."""
        path = self._path
        for i, item in enumerate(path):
            if isinstance(item, Unit):
                item._run = item._run._original
            elif callable(item):
                path[i] = item._original
    
    def debug(self):
        """Converge in debug mode. Just try it!"""
        self._debug_on()
        try: self._converge()
        finally: self._debug_off()
        end = self._error_info()
        if end:
            print(f'\nFinished debugging{end}')
        else:
            print('\n        Finished debugging')
            
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
        if bst.ALWAYS_DISPLAY_DIAGRAMS:
            try: self.diagram('minimal')
            except: pass
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

    def show(self, T=None, P=None, flow=None, composition=None, N=None, 
             IDs=None, data=True):
        """Prints information on system."""
        print(self._info(T, P, flow, composition, N, IDs, data))

    def _info(self, T, P, flow, composition, N, IDs, data):
        """Return string with all specifications."""
        error = self._error_info()
        ins_and_outs = repr_ins_and_outs(self.ins, self.outs, 
                                         T, P, flow, composition, N, IDs, data)
        return (f"System: {self.ID}"
                + error + '\n'
                + ins_and_outs)


class FacilityLoop:
    __slots__ = ('system', '_recycle',
                 'maxiter', 'molar_tolerance', 'temperature_tolerance',
                 'relative_molar_tolerance', 'relative_temperature_tolerance',
                 'alternative_convergence_check',
                 '_converge_method', '_mol_error', '_T_error', 
                 '_rmol_error', '_rT_error','_iter')
    
    default_maxiter = 50
    default_molar_tolerance = System.default_molar_tolerance
    default_temperature_tolerance = System.default_temperature_tolerance
    default_relative_molar_tolerance = System.default_relative_molar_tolerance
    default_relative_temperature_tolerance = System.default_relative_temperature_tolerance
    default_converge_method = System.default_converge_method
    
    def __init__(self, system, recycle):
        self.system = system
        self._recycle = recycle
        self._reset_errors()
        self._load_defaults()
        
    recycle = System.recycle
    converge_method = System.converge_method
    _get_recycle_temperatures = System._get_recycle_temperatures
    _get_recycle_data = System._get_recycle_data
    _set_recycle_data = System._set_recycle_data
    _reset_errors = System._reset_errors
    _error_info = System._error_info
    _iter_run = System._iter_run
    _fixedpoint = System._fixedpoint
    _aitken = System._aitken
    _wegstein = System._wegstein
    _solve = System._solve
    _load_defaults = System._load_defaults
    
    def _reset_iter(self):
        self.system._reset_iter()
        self._iter = 0
    
    def _run(self): self.system.simulate()
    
    def __call__(self): self._converge_method()

    def __repr__(self): 
        return f"<{type(self).__name__}: {self.system.ID}>"
    
       
from biosteam import _flowsheet as flowsheet_module