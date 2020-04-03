# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
from flexsolve import SolverError, conditional_wegstein, conditional_aitken
from ._digraph import (digraph_from_units_and_streams,
                       minimal_digraph,
                       surface_digraph,
                       finalize_digraph)
from thermosteam import Stream
from thermosteam.utils import registered
from ._exceptions import try_method_with_object_stamp
from ._network import Network
from ._facility import Facility
from ._unit import Unit
from ._report import save_report
from .utils import colors, strtuple
import biosteam as bst

__all__ = ('System',)

# %% Functions for building systems

def streams_from_path(path):
    isa = isinstance
    streams = set()
    for i in path:
        if isa(i, System):
            streams.add(i.streams)
        elif isa(i, Unit):
            streams.update(i._ins + i._outs)
    return streams

def feeds_from_streams(streams):
    return {s for s in streams if not s._source}

def products_from_streams(streams):
    return {s for s in streams if not s._sink}

def filter_out_missing_streams(streams):
    streams.intersection_update([i for i in streams if i])

# %% Functions for taking care of numerical specifications within a system path
    
def run_unit_in_path(unit):
    numerical_specification = unit._numerical_specification
    if numerical_specification:
        method = numerical_specification
    else:
        method = unit._run
    try_method_with_object_stamp(unit, method)

def converge_system_in_path(system):
    numerical_specification = system._numerical_specification
    if numerical_specification:
        method = numerical_specification
    else:
        method = system._converge
    try_method_with_object_stamp(system, method)

def simulate_unit_in_path(unit):
    numerical_specification = unit._numerical_specification
    if numerical_specification:
        method = numerical_specification
    else:
        method = unit.simulate
    try_method_with_object_stamp(unit, method)

def simulate_system_in_path(system):
    numerical_specification = system._numerical_specification
    if numerical_specification:
        method = numerical_specification
    else:
        method = system.simulate
    try_method_with_object_stamp(system, method)

# %% Debugging and exception handling

def _evaluate(self, command=None):
    """Evaluate a command and request user input for next command. If no command, return. This function is used for debugging a System object."""    
    # Done evaluating if no command, exit debugger if 'exit'
    if command is None:
        Next = colors.next('Next: ') + f'{repr(self)}\n'
        info = colors.info("Enter to continue or type to evaluate:\n")
        command = input(Next + info + ">>> ")
    
    if command == 'exit': raise KeyboardInterrupt()
    if command:
        # Build locals dictionary for evaluating command
        lcs = {self.ID: self, 'bst': bst} 
        try:
            out = eval(command, {}, lcs)            
        except Exception as err:
            # Print exception and ask to raise error or continue evaluating
            err = colors.exception(f'{type(err).__name__}:') + f' {str(err)}\n\n'
            info = colors.info(f"Enter to raise error or type to evaluate:\n")
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
        return self._converge.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' == method:
            self._converge = self._wegstein
        elif 'fixedpoint' == method:
            self._converge = self._fixed_point
        elif 'aitken' == method:
            self._converge = self._aitken
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
    recycle=None : :class:`~thermosteam.Stream`, optional
        A tear stream for the recycle loop.
    facilities=() : tuple[Unit, function, and/or System], optional
        Offsite facilities that are simulated only after
        completing the path simulation.

    """
    ### Class attributes ###
    
    #: Maximum number of iterations
    maxiter = 100

    #: Molar tolerance (kmol/hr)
    molar_tolerance = 1.0
    
    #: Temperature tolerance (K)
    T_tolerance = 0.10

    # [dict] Cached downstream systems by (system, unit, with_facilities) keys
    _cached_downstream_systems = {} 

    @classmethod
    def from_feedstock(cls, ID, feedstock, feeds=None, facilities=(), ends=None):
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
        
        """
        network = Network.from_feedstock(feedstock, feeds, ends)
        return cls.from_network(ID, network, facilities)

    @classmethod
    def from_network(cls, ID, network, facilities=()):
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
        
        """
        facilities = Facility.ordered_facilities(facilities)
        isa = isinstance
        path = tuple([(cls.from_network('', i) if isa(i, Network) else i)
                      for i in network.path])
        self = cls.__new__(cls)
        self.units = network.units
        self.streams = streams = network.streams
        self.feeds = feeds = network.feeds
        self.products = products = network.products
        self._numerical_specification = None
        self._reset_errors()
        self._set_path(path)
        self._set_facilities(facilities)
        self._set_recycle(network.recycle)
        self._register(ID)
        if facilities:
            f_streams = streams_from_path(facilities)
            f_feeds = feeds_from_streams(f_streams)
            f_products = products_from_streams(f_streams)
            streams.update(f_streams)
            feeds.update(f_feeds)
            products.update(f_products)
        self._finalize_streams()
        return self
        
    def __init__(self, ID, path, recycle=None, facilities=()):
        self._numerical_specification = None
        self._load_flowsheet()
        self._reset_errors()
        self._set_path(path)
        self._load_units()
        self._set_facilities(facilities)
        self._load_streams()
        self._finalize_streams()
        self._set_recycle(recycle)
        self._register(ID)
    
    numerical_specification = Unit.numerical_specification
    save_report = save_report
    
    def _load_flowsheet(self):
        self.flowsheet = flowsheet_module.main_flowsheet.get_flowsheet()
    
    def _set_recycle(self, recycle):
        if recycle is None:
            self._converge = self._run
        else:
            assert isinstance(recycle, Stream), (
             "recycle must be a Stream instance or None, not "
            f"{type(recycle).__name__}")
        self._recycle = recycle
    
    def _set_path(self, path):
        #: tuple[Unit, function and/or System] A path that is run element
        #: by element until the recycle converges.
        self.path = path
        
        #: set[System] All subsystems in the system
        self.subsystems = subsystems = set()
        
        #: list[Unit] Network of only unit operations
        self._unit_path = unit_path = []
        
        isa = isinstance
        for i in path:
            if i in unit_path: continue
            if isa(i, Unit): 
                unit_path.append(i)
            elif isa(i, System):
                unit_path.extend(i._unit_path)
                subsystems.add(i)
    
        #: set[Unit] All units in the path that have costs
        self._path_costunits = costunits = {i for i in unit_path
                                            if i._design or i._cost}
        
        #: set[Unit] All units that have costs.
        self._costunits = costunits = costunits.copy()
        
    def _load_units(self):
        #: set[Unit] All units within the system
        self.units = set(self._unit_path)
    
    def _set_facilities(self, facilities):
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the path simulation.
        self.facilities = facilities = tuple(facilities)
        subsystems = self.subsystems
        costunits = self._costunits
        units = self.units
        isa = isinstance
        for i in facilities:
            if isa(i, Unit):
                i._load_stream_links()
                units.add(i)
                if i._cost: costunits.add(i)
                if isa(i, Facility): i._system = self
            elif isa(i, System):
                units.update(i.units)
                subsystems.add(i)
                costunits.update(i._costunits)
        
    def _load_streams(self):
        #: set[:class:`~thermosteam.Stream`] All streams within the system
        self.streams = streams = set()
        
        for u in self.units:
            streams.update(u._ins + u._outs)
        for sys in self.subsystems:
            streams.update(sys.streams)
        
        #: set[:class:`~thermosteam.Stream`] All feed streams in the system.
        self.feeds = feeds_from_streams(streams)
        
        #: set[:class:`~thermosteam.Stream`] All product streams in the system.
        self.products = products_from_streams(streams)
        
    def _load_stream_links(self):
        for u in self._unit_path: u._load_stream_links()
        
    def _filter_out_missing_streams(self):
        for stream_set in (self.streams, self.feeds, self.products):
            filter_out_missing_streams(stream_set)
        
    def _finalize_streams(self):
        self._load_stream_links()
        self._filter_out_missing_streams()
        
    @property
    def TEA(self):
        """[TEA] Object for Techno-Economic Analysis."""
        try: return self._TEA
        except AttributeError: return None
    
    @property
    def recycle(self):
        """[:class:`~thermosteam.Stream`] A tear stream for the recycle loop"""
        return self._recycle

    @property
    def converge_method(self):
        """Iterative convergence method ('wegstein', 'aitken', or 'fixed point')."""
        return self._converge.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        if self.recycle is None:
            raise ValueError(
                "cannot set converge method when no recyle is specified")
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' == method:
            self._converge = self._wegstein
        elif 'fixedpoint' == method:
            self._converge = self._fixed_point
        elif 'aitken' == method:
            self._converge = self._aitken
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
                elif i in downstream_units:
                    path.append(i)
                elif (not isa(i, Unit)
                      or i.line == 'Balance'):
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
            downstream_facilities = self.facilities            
        else:
            unit_found = False
            isa = isinstance
            for pos, i in enumerate(self.facilities):
                if unit is i or (isa(i, System) and unit in i.units):
                    downstream_facilities = self.facilities[pos:]
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
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'")
        finalize_digraph(f, file=file, format=format, **graph_attrs)
            
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
        recycle = self.recycle
        rmol = recycle.mol
        rmol[:] = mol
        T = recycle.T
        self._run()
        mol_error = abs(mol - recycle.mol).sum()
        T_error = abs(T - recycle.T)
        self._iter += 1
        if mol_error < self.molar_tolerance and T_error < self.T_tolerance:
            unconverged = False
        elif self._iter > self.maxiter:
            raise SolverError(f'{repr(self)} could not converge' + self._error_info())
        else:
            unconverged = True
        self._T_error = T_error
        self._mol_error = mol_error
        return rmol.copy(), unconverged
            
        
    def _setup(self):
        """Setup each element of the system."""
        isa = isinstance
        for a in self.path:
            if isa(a, (Unit, System)): a._setup()
            else: pass # Assume it is a function
        
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
        """Converge system recycle using inner and outer loops with fixed-point iteration."""
        r = self.recycle
        rmol = r.mol
        while True:
            mol = rmol.copy()
            T = r.T
            self._run()
            self._mol_error = abs(mol - rmol).sum()
            self._T_error = abs(T - r.T)
            self._iter += 1
            if (self._mol_error < self.molar_tolerance
                and self._T_error < self.T_tolerance): break
            if self._iter > self.maxiter:
                raise SolverError(f'{repr(self)} could not converge' + self._error_info())
            
    def _wegstein(self):
        """Converge the system recycle iteratively using wegstein's method."""
        conditional_wegstein(self._iter_run, self.recycle.mol.copy())
    
    def _aitken(self):
        """Converge the system recycle iteratively using Aitken's method."""
        conditional_aitken(self._iter_run, self.recycle.mol.copy())
    
    # Default converge method
    _converge = _aitken

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
        
        #: Specification error
        self._spec_error = 0
        
        #: Number of iterations
        self._iter = 0
    
    def reset_flows(self):
        """Reset all process streams to zero flow."""
        self._reset_errors()        
        feeds = self.feeds
        for stream in self.streams:
            if stream not in feeds: stream.empty()

    def simulate(self):
        """Converge the path and simulate all units."""
        self._reset_iter()
        self._setup()
        self._converge()
        for i in self._path_costunits:
            try_method_with_object_stamp(i, i._summary)
        isa = isinstance
        for i in self.facilities:
            if isa(i, Unit):
                simulate_unit_in_path(i)
            elif isa(i, System):
                simulate_system_in_path(i)
            else:
                i() # Assume it is a function
        
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
            print(f'\n        Finished debugging')

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
    
    def _ipython_display_(self):
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
            recycle = f"\n recycle: {self.recycle}"
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


from biosteam import _flowsheet as flowsheet_module