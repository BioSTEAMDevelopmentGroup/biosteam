# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
from .digraph import (digraph_from_units_and_streams,
                      minimal_digraph,
                      surface_digraph,
                      finalize_digraph)
from thermosteam import Stream, MultiStream
from thermosteam.utils import registered
from .exceptions import try_method_with_object_stamp
from ._network import Network
from ._facility import Facility
from ._unit import Unit
from .report import save_report
from .exceptions import InfeasibleRegion
from .utils import colors, strtuple
from collections.abc import Iterable
import biosteam as bst
import numpy as np

__all__ = ('System',)    

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
    specification = unit._specification
    if specification:
        try_method_with_object_stamp(unit, unit._load_stream_links)
        try_method_with_object_stamp(unit, unit._setup)
        try_method_with_object_stamp(unit, specification)
        try_method_with_object_stamp(unit, unit._summary)
    else:
        try_method_with_object_stamp(unit, unit.simulate)

def simulate_system_in_path(system):
    specification = system._specification
    if specification:
        method = specification
    else:
        method = system.simulate
    try_method_with_object_stamp(system, method)

# %% Debugging and exception handling

def raise_recycle_type_runtime_error(recycle):
    raise RuntimeError(
       f"invalid recycle of type '{type(recycle).__name__}' encountered; "
        "recycle must be either a Stream object, a tuple of Stream objects, or None"
    )

def _evaluate(self, command=None):
    """
    Evaluate a command and request user input for next command.
    If no command, return. This function is used for debugging a System object.
    """    
    # Done evaluating if no command, exit debugger if 'exit'
    if command is None:
        Next = colors.next('Next: ') + f'{repr(self)}\n'
        info = colors.info("Enter to continue or type to evaluate:\n")
        command = input(Next + info + ">>> ")
    
    if command == 'exit': raise KeyboardInterrupt()
    if command:
        # Build locals dictionary for evaluating command
        F = bst.main_flowsheet
        lcs = {self.ID: self, 'bst': bst,
               **F.system.__dict__,
               **F.stream.__dict__,
               **F.unit.__dict__,
               **F.flowsheet.__dict__
        } 
        try:
            out = eval(command, {}, lcs)            
        except Exception as err:
            # Print exception and ask to raise error or continue evaluating
            err = colors.exception(f'{type(err).__name__}:') + f' {str(err)}\n\n'
            info = colors.info("Enter to raise error or type to evaluate:\n")
            command = input(err + info + ">>> ")
            if command == '': raise err
            _evaluate(self, command)        
        else:
            # If successful, continue evaluating
            if out is None: pass
            elif (not hasattr(out, '_ipython_display_')
                  or isinstance(out, type)): print(out)
            else: out._ipython_display_()
                
            command = input(">>> ")
            _evaluate(self, command)

def _method_debug(self, func):
    """Method decorator for debugging system."""
    def wrapper(*args, **kwargs):
        # Run method and ask to evaluate
        _evaluate(self)
        func(*args, **kwargs)
        
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper

def _notify_run_wrapper(self, func):
    """Decorate a System run method to notify you after each loop"""
    def wrapper(*args, **kwargs):
        if self.recycle:
            func(*args, **kwargs)
            input(f'        Finished loop #{self._iter}\n')
        else:
            func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


# %% Process flow

class system(type):
    @property
    def converge_method(self):
        """Iterative convergence method ('wegstein', 'aitken', or 'fixed point')."""
        return self._converge_method.__name__[1:]
    @converge_method.setter
    def converge_method(self, method):
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' == method:
            self._converge_method = self._wegstein
        elif 'fixedpoint' == method:
            self._converge_method = self._fixed_point
        elif 'aitken' == method:
            self._converge_method = self._aitken
        else:
            raise ValueError(f"only 'wegstein', 'aitken', and 'fixed point' methods are valid, not '{method}'")

@registered('SYS')
class System(metaclass=system):
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
    path : tuple[Unit, function and/or System]
        A path that is run element by element until the recycle converges.
    recycle=None : :class:`~thermosteam.Stream` or tuple[:class:`~thermosteam.Stream`], optional
        A tear stream for the recycle loop.
    facilities=() : tuple[Unit, function, and/or System], optional
        Offsite facilities that are simulated only after
        completing the path simulation.
    facility_recycle : :class:`~thermosteam.Stream`, optional
        Recycle stream between facilities and system path.
    converge_method : str, optional
        Name of iteration method to converge recycle stream.
    N_runs : int, optional
        Number of iterations to run system. This parameter is applicable 
        only to systems with no recycle loop.
    hx_convergence : 'rigorous' or 'estimate', optional
        If 'rigorous', recycle streams to process heat exchangers are rigorously
        converged. If 'estimate', loops with process heat exchanger are only 
        runned twice.

    """
    ### Class attributes ###
    
    #: Maximum number of iterations
    maxiter = 200
    
    #: Molar tolerance (kmol/hr)
    molar_tolerance = 0.50
    
    #: Temperature tolerance (K)
    temperature_tolerance = 0.10

    # [dict] Cached downstream systems by (system, unit, with_facilities) keys
    _cached_downstream_systems = {} 

    @classmethod
    def from_feedstock(cls, ID, feedstock, feeds=None, facilities=(), 
                       ends=None, facility_recycle=None, hx_convergence='rigorous'):
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
        hx_convergence : 'rigorous' or 'estimate', optional
            If 'rigorous', recycle streams to process heat exchangers are rigorously
            converged. If 'estimate', loops with process heat exchanger are only 
            runned twice.
        
        """
        network = Network.from_feedstock(feedstock, feeds, ends)
        return cls.from_network(ID, network, facilities, facility_recycle, hx_convergence)

    @classmethod
    def from_network(cls, ID, network, facilities=(), facility_recycle=None, hx_convergence='rigorous'):
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
        hx_convergence : 'rigorous' or 'estimate', optional
            If 'rigorous', recycle streams to process heat exchangers are rigorously
            converged. If 'estimate', loops with process heat exchanger are only 
            runned twice.
        
        """
        facilities = Facility.ordered_facilities(facilities)
        isa = isinstance
        path = tuple([(cls.from_network('', i) 
                       if isa(i, Network) else i) for i in network.path])
        self = cls.__new__(cls)
        self.units = network.units
        self.streams = streams = network.streams
        self.feeds = feeds = network.feeds
        self.products = products = network.products
        self._N_runs = self._specification = None
        self._set_recycle(network.recycle, hx_convergence)
        self._reset_errors()
        self._set_path(path)
        self._set_facilities(facilities)
        self._set_facility_recycle(facility_recycle)
        self._register(ID)
        if facilities:
            f_streams = bst.utils.streams_from_path(facilities)
            f_feeds = bst.utils.feeds(f_streams)
            f_products = bst.utils.products(f_streams)
            streams.update(f_streams)
            feeds.update(f_feeds)
            products.update(f_products)
        self._finalize_streams()
        return self
        
    def __init__(self, ID, path, recycle=None, facilities=(), 
                 facility_recycle=None, converge_method=None,
                 N_runs=None, hx_convergence='rigorous'):
        self._N_runs = self._specification = None
        self._set_recycle(recycle, hx_convergence)
        self._load_flowsheet()
        self._reset_errors()
        self._set_path(path)
        self._load_units()
        self._set_facilities(facilities)
        self._set_facility_recycle(facility_recycle)
        self._load_streams()
        self._finalize_streams()
        self._register(ID)
        if N_runs: self.N_runs = N_runs
        if converge_method: self.converge_method = converge_method
    
    specification = Unit.specification
    save_report = save_report
    
    def _extend_flattend_path_and_recycles(self, path, recycles, subsystems):
        isa = isinstance
        recycle = self.recycle
        if recycle:
            if isa(recycle, Stream):
                recycles.append(recycle)
            elif isa(recycle, Iterable):
                recycles.extend(recycle)
            else:
                raise_recycle_type_runtime_error(recycle)
        for i in self.path:
            if isa(i, System):
                if i.facilities:
                    warn(RuntimeWarning('subsystem with facilities could not be flattened'))
                    subsystems.add(i); path.append(i)
                elif i.specification:
                    warn(RuntimeWarning('subsystem with specification could not be flattened'))
                    subsystems.add(i); path.append(i)
                else:
                    i._extend_flattend_path_and_recycles(path, recycles, subsystems)
            else:
                path.append(i)
    
    def flatten(self, hx_convergence='rigorous'):
        """
        Flatten system by removing subsystems.

        Parameters
        ----------
        hx_convergence : 'rigorous' or 'estimate', optional
            If 'rigorous', recycle streams to process heat exchangers are rigorously
            converged. If 'estimate', loops with process heat exchanger are only 
            runned twice.

        """
        recycles = []
        path = []
        subsystems = set()
        self._extend_flattend_path_and_recycles(path, recycles, subsystems)
        self.path = tuple(path)
        self._set_recycle(recycles, hx_convergence)
        N_recycles = len(recycles)
        self.molar_tolerance *= N_recycles
        self.temperature_tolerance *= N_recycles
    
    @property
    def N_runs(self):
        return self._N_runs
    @N_runs.setter
    def N_runs(self, N_runs):
        if self._recycle:
            raise AttributeError("cannot set number of runs when a recycle is defined")
        self._N_runs = N_runs
    
    def _load_flowsheet(self):
        self.flowsheet = flowsheet_module.main_flowsheet.get_flowsheet()
    
    def _set_recycle(self, recycle, hx_convergence):
        isa = isinstance
        if recycle is None:
            self._recycle = recycle
        elif isa(recycle, Stream):
            if isa(recycle.sink, bst.HXprocess):
                if hx_convergence == 'rigorous':
                    self._recycle = recycle
                elif hx_convergence == 'estimate':
                    self._N_runs = 2
                    self._recycle = None
                else:
                    raise ValueError("hx_convergence must be either 'rigorous' "
                                     "or 'estimate'")
            else:
                self._recycle = recycle
        elif isa(recycle, Iterable) and all([isa(i, Stream) for i in recycle]):
            recycle = tuple(recycle)
            if all([isa(i.sink, bst.HXprocess) for i in recycle]):
                if hx_convergence == 'rigorous':
                    self._recycle = recycle
                elif hx_convergence == 'estimate':
                    self._N_runs = 2
                    self._recycle = None
                else:
                    raise ValueError("hx_convergence must be either 'rigorous' "
                                     "or 'estimate'")
            else:
                self._recycle = recycle
        else:
            raise ValueError("recycle must be either a Stream object, a tuple of Stream objects, or None; "
                            f"not {type(recycle).__name__}")
    
    def _set_path(self, path):
        #: tuple[Unit, function and/or System] A path that is run element
        #: by element until the recycle converges.
        self.path = path
        
        #: set[System] All subsystems in the system
        self.subsystems = subsystems = set()
        
        #: list[Unit] Network of only unit operations
        self._unit_path = unit_path = []
        
        #: set[Unit] All units that have costs.
        self._costunits = costunits = set()
        
        isa = isinstance
        for i in path:
            if i in unit_path: continue
            if isa(i, Unit): 
                unit_path.append(i)
            elif isa(i, System):
                unit_path.extend(i._unit_path)
                subsystems.add(i)
                costunits.update(i._costunits)
    
        #: set[Unit] All units in the path that have costs
        self._path_costunits = path_costunits = {i for i in unit_path
                                                 if i._design or i._cost}
        costunits.update(path_costunits)
        
        
    def _load_units(self):
        #: set[Unit] All units within the system
        self.units = set(self._unit_path) | self._costunits
    
    def _set_facilities(self, facilities):
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the path simulation.
        self._facilities = facilities = tuple(facilities)
        subsystems = self.subsystems
        costunits = self._costunits
        units = self.units
        isa = isinstance
        for i in facilities:
            if isa(i, Unit):
                i._load_stream_links()
                units.add(i)
                if i._cost: costunits.add(i)
                if isa(i, Facility) and not i._system: i._system = self
            elif isa(i, System):
                units.update(i.units)
                subsystems.add(i)
                costunits.update(i._costunits)
    
    def _set_facility_recycle(self, recycle):
        if recycle:
            system = self._downstream_system(recycle.sink)
            #: [FacilityLoop] Recycle loop for converging facilities
            self._facility_loop = FacilityLoop(system, recycle)
        else:
            self._facility_loop = None
        
    def _load_streams(self):
        #: set[:class:`~thermosteam.Stream`] All streams within the system
        self.streams = streams = set()
        
        for u in self.units:
            streams.update(u._ins + u._outs)
        for sys in self.subsystems:
            streams.update(sys.streams)
        
        #: set[:class:`~thermosteam.Stream`] All feed streams in the system.
        self.feeds = bst.utils.feeds(streams)
        
        #: set[:class:`~thermosteam.Stream`] All product streams in the system.
        self.products = bst.utils.products(streams)
        
    def _load_stream_links(self):
        for u in self._unit_path: u._load_stream_links()
        
    def _filter_out_missing_streams(self):
        for stream_set in (self.streams, self.feeds, self.products):
            bst.utils.filter_out_missing_streams(stream_set)
        
    def _finalize_streams(self):
        self._load_stream_links()
        self._filter_out_missing_streams()
        
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
        """[:class:`~thermosteam.Stream`] A tear stream for the recycle loop"""
        return self._recycle

    @property
    def converge_method(self):
        """Iterative convergence method ('wegstein', 'aitken', or 'fixed point')."""
        return self._converge_method.__name__[1:]
    @converge_method.setter
    def converge_method(self, method):
        if self.recycle is None:
            raise ValueError(
                "cannot set converge method when no recyle is specified")
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' == method:
            self._converge_method = self._wegstein
        elif 'fixedpoint' == method:
            self._converge_method = self._fixed_point
        elif 'aitken' == method:
            self._converge_method = self._aitken
        else:
            raise ValueError(
                f"only 'wegstein', 'aitken', and 'fixed point' methods "
                f"are valid, not '{method}'")

    
    def _downstream_path(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.units: return []
        elif self._recycle: return self.path
        unit_found = False
        downstream_units = unit._downstream_units
        path = []
        isa = isinstance
        for i in self.path:
            if unit_found:
                if isa(i, System):
                    for u in i.units:
                        if u in downstream_units:
                            path.append(i)
                            break
                elif i in downstream_units or not isa(i, Unit):
                    path.append(i)
            else:
                if unit is i:
                    unit_found = True
                    path.append(unit)
                elif isa(i, System) and unit in i.units:
                    unit_found = True   
                    path.append(i)
        return path

    def _downstream_system(self, unit):
        """Return a system with a path composed of the `unit` and
        everything downstream (facilities included)."""
        if unit is self.path[0]: return self
        system = self._cached_downstream_systems.get((self, unit))
        if system: return system
        path = self._downstream_path(unit)
        if path:
            downstream_facilities = self._facilities            
        else:
            unit_found = False
            isa = isinstance
            for pos, i in enumerate(self._facilities):
                if unit is i or (isa(i, System) and unit in i.units):
                    downstream_facilities = self._facilities[pos:]
                    unit_found = True
                    break
            assert unit_found, f'{unit} not found in system'
        system = System(None, path,
                        facilities=downstream_facilities)
        system._ID = f'{type(unit).__name__}-{unit} and downstream'
        self._cached_downstream_systems[unit] = system
        return system
    
    def _minimal_digraph(self, **graph_attrs):
        """Return digraph of the path as a box."""
        return minimal_digraph(self.ID, self.units, self.streams, **graph_attrs)

    def _surface_digraph(self, **graph_attrs):
        return surface_digraph(self.path)

    def _thorough_digraph(self, **graph_attrs):
        return digraph_from_units_and_streams(self.units, self.streams, 
                                              **graph_attrs)
        
    def diagram(self, kind='surface', file=None, format='png', **graph_attrs):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the system.
        
        Parameters
        ----------
        kind='surface' : {'thorough', 'surface', 'minimal'}:
            * **'thorough':** Display every unit within the path.
            * **'surface':** Display only elements listed in the path.
            * **'minimal':** Display path as a box.
        file=None : str, display in console by default
            File name to save diagram.
        format='png' : str
            File format (e.g. "png", "svg").
        
        """
        if kind == 'thorough':
            f = self._thorough_digraph(format=format, **graph_attrs)
        elif kind == 'surface':
            f = self._surface_digraph(format=format, **graph_attrs)
        elif kind == 'minimal':
            f = self._minimal_digraph(format=format, **graph_attrs)
        else:
            raise ValueError("kind must be either 'thorough', 'surface', or 'minimal'")
        finalize_digraph(f, file, format)
            
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
        rmol : numpy.ndarray
               New recycle molar flow rates.
        unconverged : bool
                      True if recycle has not converged.
            
        """
        if (mol < 0.).any():
            raise InfeasibleRegion('material flow')
        self._set_recycle_data(mol)
        T = self._get_recycle_temperatures()
        self._run()
        rmol = self._get_recycle_data()
        rT = self._get_recycle_temperatures()
        self._mol_error = mol_error = abs(mol - rmol).sum()
        self._T_error = T_error = np.sum(abs(T - rT))
        self._iter += 1
        if mol_error < self.molar_tolerance and T_error < self.temperature_tolerance:
            unconverged = False
        elif self._iter == self.maxiter:
            raise RuntimeError(f'{repr(self)} could not converge' + self._error_info())
        else:
            unconverged = True
        return rmol, unconverged
            
    def _get_recycle_data(self):
        recycle = self.recycle
        if isinstance(recycle, Stream):
            return recycle.imol.data.copy()
        elif isinstance(recycle, Iterable):
            return np.vstack([i.imol.data for i in recycle])
        else:
            raise RuntimeError('no recycle available')
    
    def _set_recycle_data(self, data):
        recycle = self.recycle
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
        recycle = self.recycle
        if isinstance(recycle, Stream):
            return self.recycle.T
        elif isinstance(recycle, Iterable):
            return np.array([i.T for i in recycle], float)
        else:
            raise RuntimeError('no recycle available')
    
    def _setup(self):
        """Setup each element of the system."""
        isa = isinstance
        for i in self.path:
            if isa(i, (Unit, System)): i._setup()
        
    def _run(self):
        """Rigorous run each element of the system."""
        isa = isinstance
        for i in self.path:
            if isa(i, Unit):
                run_unit_in_path(i)
            elif isa(i, System): 
                converge_system_in_path(i)
            else: i() # Assume it is a function
    
    # Methods for convering the recycle stream    
    def _fixed_point(self):
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
    
    # Default converge method
    _converge_method = _aitken
   
    def _converge(self):
        if self._N_runs:
            for i in range(self.N_runs): self._run()
        elif self._recycle:
            self._converge_method()
        else:
            self._run()
        
    def _design_and_cost(self):
        for i in self._path_costunits:
            try_method_with_object_stamp(i, i._summary)
        isa = isinstance
        for i in self._facilities:
            if isa(i, Unit):
                simulate_unit_in_path(i)
            elif isa(i, System):
                simulate_system_in_path(i)
            else:
                i() # Assume it is a function

    def _reset_iter(self):
        self._iter = 0
        for system in self.subsystems: system._reset_iter()
    
    def reset_names(self, unit_format=None, stream_format=None):
        """Reset names of all streams and units according to the path order."""
        Unit._default_ID = unit_format if unit_format else ['U', 0]
        Stream._default_ID = stream_format if stream_format else ['d', 0]
        streams = set()
        units = set()
        for i in self._unit_path:
            if i in units: continue
            try: i.ID = ''
            except: continue
            for s in (i._ins + i._outs):
                if (s and s._sink and s._source
                    and s not in streams):
                    s.ID = ''
                    streams.add(s)
            units.add(i)
    
    def _reset_errors(self):
        #: Molar flow rate error (kmol/hr)
        self._mol_error = 0
        
        #: Temperature error (K)
        self._T_error = 0
        
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
        recycle = self.recycle
        if recycle:
            if isinstance(recycle, Stream):
                recycle.empty()
            elif isinstance(recycle, Iterable):
                for i in recycle: i.empty()
            else:
                raise_recycle_type_runtime_error(recycle)
        for system in self.subsystems:
            system.empty_recycles()

    def reset_cache(self):
        """Reset cache of all unit operations."""
        for unit in self.units: unit.reset_cache()

    def simulate(self):
        """Converge the path and simulate all units."""
        self._setup()
        self._converge()
        self._design_and_cost()
        if self._facility_loop: self._facility_loop()
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        self._run = _notify_run_wrapper(self, self._run)
        self.path = path = list(self.path)
        for i, item in enumerate(path):
            if isinstance(item, Unit):
                item._run = _method_debug(item, item._run)
            elif isinstance(item, System):
                item._converge = _method_debug(item, item._converge)
            elif callable(item):
                path[i] = _method_debug(item, item)

    def _debug_off(self):
        """Turn off debug mode."""
        self._run = self._run._original
        path = self.path
        for i, item in enumerate(path):
            if isinstance(item, Unit):
                item._run = item._run._original
            elif isinstance(item, System):
                item._converge = item._converge._original
            elif callable(item):
                path[i] = item._original
        self.path = tuple(path)
    
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
    def __str__(self):
        if self.ID: return self.ID
        else: return type(self).__name__ 
    
    def __repr__(self):
        if self.ID: return f'<{type(self).__name__}: {self.ID}>'
        else: return f'<{type(self).__name__}>'

    def show(self):
        """Prints information on unit."""
        print(self._info())
    
    def to_network(self):
        """Return network that defines the system path."""
        isa = isinstance
        path = [(i.to_network() if isa(i, System) else i) for i in self.path]
        network = Network.__new__(Network)    
        network.path = path
        network.recycle = self.recycle
        network.units = self.units
        network.subnetworks = [i for i in path if isa(i, Network)]
        network.feeds = self.feeds
        network.products = self.products
        return network
    
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
        path_info = get_path_info(self.path)
        info = update_info(path_info)
        facilities = self.facilities
        if facilities:
            facilities_info = get_path_info(self.facilities)
            facilities_info = f'facilities={facilities_info}'
            info = update_info(facilities_info)
        recycle = self.recycle
        if recycle:
            recycle = self._get_recycle_info()
            info = update_info(f"recycle={recycle}")
        if self.N_runs:
            info = update_info(f"N_runs={self.N_runs}")
        info += ')'
        return info
    
    def _get_recycle_info(self):
        recycle = self.recycle
        if isinstance(recycle, Stream):
            recycle = recycle._source_info()
        else:
            recycle = ", ".join([i._source_info() for i in recycle])
            recycle = '[' + recycle + ']'
        return recycle
    
    def _ipython_display_(self):
        if bst.ALWAYS_DISPLAY_DIAGRAMS:
            try: self.diagram('minimal')
            except: pass
        self.show()

    def _error_info(self):
        """Return information on convergence."""
        if self.recycle:
            return (f"\n convergence error: Flow rate   {self._mol_error:.2e} kmol/hr"
                    f"\n                    Temperature {self._T_error:.2e} K"
                    f"\n iterations: {self._iter}")
        else:
            return ""

    def _info(self):
        """Return string with all specifications."""
        if self.recycle is None:
            recycle = ''
        else:
            recycle = f"\n recycle: {self._get_recycle_info()}"
        error = self._error_info()
        path = strtuple(self.path)
        i = 1; last_i = 0
        while True:
            i += 2
            i = path.find(', ', i)
            i_next = path.find(', ', i+2)
            if (i_next-last_i) > 35:
                path = (path[:i] + '%' + path[i:])
                last_i = i
            elif i == -1: break
        path = path.replace('%, ', ',\n' + ' '*8)
        
        if self.facilities:
            facilities = strtuple(self.facilities)
            i = 1; last_i = 0
            while True:
                i += 2
                i = facilities.find(', ', i)
                if (i - last_i) > 35:
                    facilities = (facilities[:i] + '%' + facilities[i:])
                    last_i = i
                elif i == -1: break
            facilities = facilities.replace('%, ', ',\n'+' '*14)
            facilities = f"\n facilities: {facilities}" 
        else:
            facilities = ''
        
        return (f"System: {self.ID}"
                + recycle
                + f"\n path: {path}"
                + facilities
                + error)


class FacilityLoop(metaclass=system):
    __slots__ = ('system', 'recycle',
                 '_mol_error', '_T_error', '_iter')
    
    #: Maximum number of iterations to solve facilities
    maxiter = 50
    
    #: Molar tolerance (kmol/hr)
    molar_tolerance = 0.50
    
    #: Temperature tolerance (K)
    temperature_tolerance = 0.10
    
    def __init__(self, system, recycle):
        self.system = system
        self.recycle = recycle
        self._reset_errors()
        
    _get_recycle_temperatures = System._get_recycle_temperatures
    _get_recycle_data = System._get_recycle_data
    _set_recycle_data = System._set_recycle_data
    _reset_errors = System._reset_errors
    _error_info = System._error_info
    _iter_run = System._iter_run
    _fixed_point = System._fixed_point
    _aitken = System._aitken
    _wegstein = System._wegstein
    _converge_method = System._converge_method
    _solve = System._solve
    converge_method = System.converge_method
    
    def _reset_iter(self):
        self.system._reset_iter()
        self._iter = 0
    
    def _run(self): self.system.simulate()
    def __call__(self): self._converge_method()

    def __repr__(self): 
        return f"<{type(self).__name__}: {self.system.ID}>"
    
        
from biosteam import _flowsheet as flowsheet_module