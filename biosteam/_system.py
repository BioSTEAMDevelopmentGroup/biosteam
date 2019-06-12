# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
from copy import copy
import numpy as np
from scipy.optimize import newton
from ._exceptions import Stop, SolverError, _try_method
from ._flowsheet import find, make_digraph, save_digraph
from ._stream import Stream
from ._unit import Unit
from ._report import save_report
from ._utils import color_scheme, missing_stream, strtuple, function
from warnings import warn
import biosteam

__all__ = ('System',)

# %% Debugging and exception handling

def _evaluate(self, command=None):
    """Evaluate a command and request user input for next command. If no command, return. This function is used for debugging a System object."""    
    CS = color_scheme
    # Done evaluating if no command, exit debugger if 'exit'
    if command is None:
        Next = CS.next('Next: ') + f'{repr(self)}\n'
        info = CS.info("Enter to continue or type to evaluate:\n")
        command = input(Next + info + ">>> ")
    
    if command == 'exit': raise Stop
    if command:
        # Build locals dictionary for evaluating command
        lcs = {} 
        for attr in ('stream', 'unit', 'system'):
            dct = getattr(find, attr).__dict__
            lcs.update({i:j() for i, j in dct.items()})
        lcs.update({attr:getattr(biosteam, attr) for attr in biosteam.__all__})
        try:
            out = eval(command, {}, lcs)            
        except Exception as err:
            # Print exception and ask to raise error or continue evaluating
            err = CS.exception(f'{type(err).__name__}:') + f' {str(err)}\n\n'
            info = CS.info(f"Enter to raise error or type to evaluate:\n")
            command = input(err + info + ">>> ")
            if command == '': raise err
            _evaluate(self, command)        
        else:
            # If successful, continue evaluating
            if out is None: pass
            elif (not hasattr(out, '_ipython_display_')
                  or isinstance(type(out), type)): print(out)
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
            x = self._solver_error['iter']
            input(f'        Finished loop #{x}\n')
        else:
            func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper
    
# %% System node for diagram

class _systemUnit(Unit):
    """Dummy unit for displaying a system."""
    line = 'System'
    ID = None
    _ID = property(lambda self: self.ID)

_sysgraphics = _systemUnit._graphics
_sysgraphics.edge_in = _sysgraphics.edge_in * 10
_sysgraphics.edge_out = _sysgraphics.edge_out * 15
_sysgraphics.node['peripheries'] = '2'

class _streamUnit(Unit):
    line = ''
    ID = None
    _ID = _systemUnit._ID
    
_stream_graphics = _streamUnit._graphics
_stream_graphics.node['fillcolor'] = 'white:#79dae8'
del _stream_graphics, _sysgraphics

# %% Other

_isfeed = lambda stream: not stream._source and stream._sink
_isproduct = lambda stream: not stream._sink and stream._source


# %% Process flow

class System:
    """Create a System object that can iteratively run each element in a network of BioSTREAM objects until the recycle stream is converged. A network can have function, Unit and/or System objects. When the network contains an inner System object, it converges/solves it in each loop/iteration.

    **Parameters**

         **ID:** [str] A unique identification. If ID is None, instance will not be registered in flowsheet.

         **network:** tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.

         **recycle:** [Stream] A tear stream for the recycle loop.
         
         **facilities:** tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.

    """
    ### Class attributes ###
    
    #: [dict] Dictionary of convergence options regarding maximum number of iterations and molar flow rate and temperature tolerances
    options = {'Maximum iteration': 100,
               'Molar tolerance (kmol/hr)': 0.01,
               'Temperature tolerance (K)': 0.10}

    # [dict] Cached downstream systems by (system, unit, with_facilities) keys
    _cached_downstream_systems = {} 

    def __init__(self, ID, network, recycle=None, facilities=None):
        #: [dict] Current molar flow and temperature errors and number of iterations made
        self._solver_error = {'mol_error': 0,
                              'T_error': 0,
                              'spec_error': 0,
                              'iter': 0}
         
        #: set[Stream] All streams within the system
        self.streams = streams = set()

        #: set[System] All subsystems in the system
        self.subsystems = subsystems = set()
        
        #: list[Unit] Network of only unit operations
        self._unitnetwork = units = []
        inst = isinstance
        for i in network:
            if i in units: continue
            if inst(i, Unit): 
                units.append(i)
                streams.update(i._ins)
                streams.update(i._outs)
            elif inst(i, System):
                units.extend(i._unitnetwork)
                subsystems.add(i)
                streams.update(i.streams)
        streams.discard(missing_stream) 
        
        #: tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.
        self.network = tuple(network)
        
        # link all unit operations with linked streams
        for u in units:
            try:
                if u._linkedstreams: u._link_streams() 
            except Exception as Error:
                if missing_stream in (u._ins + u._outs):
                    warn(f'missing stream object in {repr(u)}')
        
        #: set[Unit] All units within the system
        self.units = units = set(units)
        
        #: set[Unit] All units in the network that have costs
        self._network_costunits = costunits = {i for i in units if i._has_cost}
        
        #: set[Unit] All units that have costs.
        self._costunits = costunits = costunits.copy()
        if facilities:
            for i in facilities:
                if inst(i, Unit):
                    units.add(i)
                    streams.update(i._ins + i._outs)
                    if i._has_cost: costunits.add(i)
                elif inst(i, System):
                    units.update(i.units)
                    streams.update(i.streams)
                    subsystems.add(i)
                    costunits.update(i._costunits)
            #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.
            self.facilities = tuple(facilities)
        else: self.facilities = ()
        
        has = hasattr
        upstream_connections = set()
        for s in streams:
            if has(s, '_downstream_connection'):
                upstream_connections.add(s)
        streams.difference_update(upstream_connections)
        
        #: set[Stream] All feed streams in the system.
        self.feeds = set(filter(_isfeed, streams))
        
        #: set[Stream] All product streams in the system.
        self.products = set(filter(_isproduct, streams)) 
        
        #: [TEA] System object for Techno-Economic Analysis.
        self._TEA = None
        self.recycle = recycle
        
        if ID:
            ID = ID.replace(' ', '_')
            ID_words = ID.split('_')
            if not all(word.isalnum() for word in ID_words):
                raise ValueError('ID cannot have any special characters')
            self._ID = ID
            find.system[ID] = self
    
    save_report = save_report
    
    @property
    def ID(self):
        """Identification."""
        return self._ID
    
    @property
    def TEA(self):
        """[TEA] System object for Techno-Economic Analysis."""
        return self._TEA
    
    @property
    def recycle(self):
        """A tear stream for the recycle loop"""
        return self._recycle

    @recycle.setter
    def recycle(self, stream):
        if stream is None:
            self._converge_method = self._run
        elif not isinstance(stream, Stream):
            raise ValueError(f"recycle must be a Stream instance or None, not {type(stream).__name__}")
        self._recycle = stream

    @property
    def converge_method(self):
        """Iterative convergence method ('Wegstein', or 'fixed point')."""
        return self._converge_method.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        if self.recycle is None:
            raise ValueError("cannot set converge method when no recyle is specified")
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' in method:
            self._converge_method = self._Wegstein
        elif 'fixedpoint' in method:
            self._converge_method = self._fixed_point
        else:
            raise ValueError(f"only 'Wegstein' and 'fixed point' methods are valid, not '{method}'")

    
    def _downstream_network(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.units: return []
        elif self._recycle: return self.network
        unit_found = False
        downstream_units = unit._downstream_units
        network = []
        inst = isinstance
        for i in self.network:
            if unit_found:
                if inst(i, System):
                    for u in i.units:
                        if u in downstream_units:
                            network.append(i)
                            break
                elif i in downstream_units:
                    network.append(i)
                elif (not inst(i, Unit)
                      or i.line == 'Balance'):
                    network.append(i)
            else:
                if unit is i:
                    unit_found = True
                    network.append(unit)
                elif inst(i, System) and unit in i.units:
                        unit_found = True   
                        network.append(i)
        return network

    def _downstream_system(self, unit):
        """Return a system with a network composed of the `unit` and everything downstream (facilities included)."""
        if unit is self.network[0]: return self
        system = self._cached_downstream_systems.get((self, unit))
        if system: return system
        network = self._downstream_network(unit)
        if network:
            downstream_facilities = self.facilities            
        else:
            unit_not_found = True
            inst = isinstance
            for pos, i in enumerate(self.facilities):
                if unit is i or (inst(i, System) and unit in i.units):
                    downstream_facilities = self.facilities[pos:]
                    unit_not_found = False
                    break
            if unit_not_found: raise ValueError(f'{unit} not found in system')
        system = System(None, network,
                        facilities=downstream_facilities)
        system._ID = f'{type(unit).__name__}-{unit} and downstream'
        self._cached_downstream_systems[unit] = system
        return system
    
    def _minimal_diagram(self, file, format='svg'):
        """Minimally display the network as a box."""
        outs = []
        ins = []
        for s in self.streams:
            source = s._source
            sink = s._sink
            if source in self.units and sink not in self.units:
                outs.append(s)
            elif sink in self.units and source not in self.units:
                ins.append(s)
        product = Stream(None)
        product._ID = ''
        feed = Stream(None)
        feed._ID = ''
        _streamUnit('\n'.join([i.ID for i in ins]),
                    feed)
        _streamUnit('\n'.join([i.ID for i in outs]),
                    None, product)
        unit = _systemUnit(self.ID, product, feed)
        unit.diagram(1, file, format)

    def _surface_diagram(self, file, format='svg'):
        """Display only surface elements listed in the network."""
        # Get surface items to make nodes and edges
        units = set()  
        refresh_units = set()
        for i in self.network:
            if isinstance(i, Unit):
                units.add(i)
            elif isinstance(i, System):
                outs = []
                ins = []
                feeds = []
                products = []
                for s in i.streams:
                    source = s._source
                    sink = s._sink
                    if source in i.units and sink not in i.units:
                        if sink: outs.append(s)
                        else: products.append(s)
                        refresh_units.add(source)
                    elif sink in i.units and source not in i.units:
                        if source: ins.append(s)
                        else: feeds.append(s)
                        refresh_units.add(sink)
                
                if len(feeds) > 1:
                    feed = Stream(None)
                    feed._ID = ''
                    units.add(_streamUnit('\n'.join([i.ID for i in feeds]), feed))
                    ins.append(feed)
                else: ins += feeds
                
                if len(products) > 1:
                    product = Stream(None)
                    product._ID = ''
                    units.add(_streamUnit('\n'.join([i.ID for i in products]),
                                          None, product))
                    outs.append(product)
                else: outs += products
                
                subsystem_unit = _systemUnit(i.ID, outs, ins)
                units.add(subsystem_unit)
                
        System(None, units)._thorough_diagram(file, format)
        # Reconnect how it was
        for u in refresh_units:
            u._ins[:] = u._ins
            u._outs[:] = u._outs
      
    def _thorough_diagram(self, file=None, format='svg'):
        """Thoroughly display every unit within the network."""
        # Create a digraph and set direction left to right
        f = make_digraph(self.units, self.streams)
        save_digraph(f, file, format)
        
    def diagram(self, kind='surface', file=None, format='svg'):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the system.
        
        **Parameters**
        
            **kind:** Must be one of the following:
                * **'thorough':** Thoroughly display every unit within the network
                * **'surface':** Display only surface elements listed in the network
                * **'minimal':** Minimally display the network as a box
        
            **file:** Must be one of the following:
                * [str] File name to save diagram.
                * [None] Display diagram in console.
        
            **format:** File format.
        
        """
        if kind == 'thorough':
            return self._thorough_diagram(file, format)
        elif kind == 'surface':
            return self._surface_diagram(file, format)
        elif kind == 'minimal':
            return self._minimal_diagram(file, format)
        else:
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'")
            

    # Methods for running one iteration of a loop
    def _run(self):
        """Rigorous run each element of the system."""
        inst = isinstance
        _try = _try_method
        for a in self.network:
            if inst(a, Unit): _try(a._run)
            elif inst(a, System): a._converge()
            else: a() # Assume it is a function
        self._solver_error['iter'] += 1
    
    # Methods for convering the recycle stream
    def _fixed_point(self):
        """Converge system recycle using inner and outer loops with fixed-point iteration."""
        # Reused attributes
        recycle = self.recycle
        run = self._run
        solver_error = self._solver_error
        maxiter, mol_tol, T_tol = self.options.values()

        # Outer loop
        abs_ = abs
        while True:
            mol_old = copy(recycle.mol)
            T_old = recycle.T
            run()
            mol_error = abs_(recycle.mol - mol_old).sum()
            T_error = abs_(recycle.T - T_old)
            if T_error < T_tol and mol_error < mol_tol: break
            if solver_error['iter'] > maxiter:
                solver_error['mol_error'] = mol_error
                solver_error['T_error'] = T_error
                raise SolverError(f'could not converge'
                                  + self._error_info())

        solver_error['mol_error'] = mol_error
        solver_error['T_error'] = T_error

    def _Wegstein(self):
        """Converge the system recycle iteratively using Wegstein's method."""
        # Reused attributes
        recycle = self.recycle
        run = self._run
        maxiter, mol_tol, T_tol = self.options.values()
        solver_error = self._solver_error

        # Prepare variables
        len_ = recycle._species._Nspecies + 1
        x0 = np.zeros(len_)
        gx0 = np.zeros(len_)
        x1 = np.zeros(len_)
        gx1 = np.zeros(len_)
        ones = np.ones(len_)
        s = np.zeros(len_)

        # First run
        x0[:-1] = recycle.mol
        x0[-1] = recycle.T
        run()
        x1[:-1] = gx0[:-1] = recycle.mol
        x1[-1] = gx0[-1] = recycle.T

        # Check convergence
        abs_ = abs
        mol_error = abs_(gx0[:-1] - x0[:-1]).sum()
        T_error = abs_(gx0[-1] - x0[-1])
        converged = mol_error < mol_tol and T_error < T_tol

        # Outer loop
        while not converged:
            run()
            gx1[:-1] = recycle.mol
            gx1[-1] = recycle.T

            # Check if converged
            mol_error = abs_(gx1[:-1] - x1[:-1]).sum()
            T_error = abs_(gx1[-1] - x1[-1])
            converged = mol_error < mol_tol and T_error < T_tol

            if solver_error['iter'] > maxiter:
                solver_error['mol_error'] = mol_error
                solver_error['T_error'] = T_error
                raise SolverError(f'{self} could not converge'
                                  + self._error_info())

            # Get relaxation factor and set next iteration
            x_diff = x1 - x0
            pos = x_diff != 0
            s[pos] = (gx1[pos] - gx0[pos])/x_diff[pos]
            x0 = x1
            gx0 = gx1
            s[s > 0.9] = 0.9
            w = ones/(ones-s)
            x1 = w*gx1 + (1-w)*x1
            recycle._mol[:] = x1[:-1]
            recycle.T = x1[-1]
        solver_error['mol_error'] = mol_error
        solver_error['T_error'] = T_error
    
    # Default converge method
    _converge_method = _Wegstein
    
    def _converge(self):
        """Converge the system recycle using an iterative solver (Wegstein by default)."""
        try:
            return self._converge_method()
        except Exception as e:
            if self._recycle: self._recycle.empty()
            raise e

    def set_spec(self, getter, setter, solver=newton, **kwargs):
        """Wrap a solver around the converge method.

        **Parameters**

             **getter:** [function] Returns objective value.

             **setter:** [function] Changes independent variable.
             
             **solver:** [function] Solves objective function.

             **kwargs:** [dict] Key word arguments passed to solver.

        """
        converge = self._converge_method
        solver_error = self._solver_error

        def error(val):
            setter(val)
            converge()
            solver_error['spec_error'] = e = getter()
            return e

        def converge_spec():
            x = solver(error, **kwargs)
            if solver is newton: kwargs['x0'] = x

        self._converge_method = converge_spec

    def _reset_iter(self):
        self._solver_error['iter'] = 0
        for system in self.subsystems: system._reset_iter()
    
    def reset_names(self, unit_format=None, stream_format=None):
        """Reset names of all streams and units according to the network order."""
        Unit._default_ID = unit_format if unit_format else ['U', 0]
        Stream._default_ID = stream_format if stream_format else ['d', 0]
        streams = set()
        units = set()
        for i in self._unitnetwork:
            if i in units: continue
            try: i.ID = ''
            except: continue
            for s in (i._ins + i._outs):
                if (s is not missing_stream
                    and s._sink and s._source
                    and s not in streams):
                    s.ID = ''
                    streams.add(s)
            units.add(i)
    
    def reset_flows(self):
        """Reset all process streams to zero flow."""
        self._solver_error = {'mol_error': 0,
                              'T_error': 0,
                              'spec_error': 0,
                              'iter': 0}
        feeds = self.feeds
        for stream in self.streams:
            if stream._link: continue
            if stream not in feeds: stream.empty()

    def simulate(self):
        """Converge the network and simulate all units."""
        self._reset_iter()
        self._converge()
        _try = _try_method
        for u in self._network_costunits: _try(u._summary)
        inst = isinstance
        for i in self.facilities:
            if inst(i, (Unit, System)): i.simulate()
            else: i() # Assume it is a function
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        self._run = _notify_run_wrapper(self, self._run)
        self.network = network = list(self.network)
        for i, item in enumerate(network):
            if isinstance(item, Unit):
                item._run = _method_debug(item, item._run)
            elif isinstance(item, System):
                item._converge = _method_debug(item, item._converge)
            elif isinstance(item, function):
                network[i] = _method_debug(item, item)

    def _debug_off(self):
        """Turn off debug mode."""
        self._run = self._run._original
        network = self.network
        for i, item in enumerate(network):
            if isinstance(item, Unit):
                item._run = item._run._original
            elif isinstance(item, System):
                item._converge = item._converge._original
            elif isinstance(item, function):
                network[i] = item._original
        self.network = tuple(network)
    
    def debug(self):
        """Convege in debug mode. Just try it!"""
        self._debug_on()
        try:
            self._converge_method()
        except Stop:
            self._debug_off()
        except Exception as error:
            self._debug_off()
            raise error
        else:
            end = self._error_info()
            if end:
                print(f'\nFinished debugging{end}')
            else:
                print(f'\n        Finished debugging')
            self._debug_off()

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
    
    def _ipython_display_(self):
        try: self.diagram('minimal')
        except: pass
        self.show()

    def _error_info(self):
        """Return information on convergence."""
        x = self._solver_error
        recycle = self.recycle
        error_info = ''
        spec = x['spec_error']
        if spec:
            error_info += f'\n specification error: {spec:.3g}'
        if recycle:
            error_info += f"\n convergence error: Flow rate   {x['mol_error']:.2e} kmol/hr"
            error_info += f"\n                    Temperature {x['T_error']:.2e} K"
        if spec or recycle:
            error_info += f"\n iterations: {x['iter']}"
        return error_info

    def _info(self):
        """Return string with all specifications."""
        if self.recycle is None:
            recycle = ''
        else:
            recycle = f"\n recycle: {self.recycle}"
        error = self._error_info()
        network = strtuple(self.network)
        i = 1; last_i = 0
        while True:
            i += 2
            i = network.find(', ', i)
            i_next = network.find(', ', i+2)
            if (i_next-last_i) > 35:
                network = (network[:i] + '%' + network[i:])
                last_i = i
            elif i == -1: break
        network = network.replace('%, ', ',\n'+' '*11)
        
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
                + f"\n network: {network}"
                + facilities
                + error)
