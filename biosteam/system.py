# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
from copy import copy
import IPython
from scipy.optimize import brentq, newton
from graphviz import Digraph
from .exceptions import Stop, SolverError, notify_error
from .flowsheet import find
from .stream import Stream, mol_units, T_units
from .unit import Unit
from . import np
from .utils import color_scheme, missing_stream, strtuple, function
from bookkeep import SmartBook
from .sim import Block, Grid
CS = color_scheme

__all__ = ('System',)

# %% Debugging and exception handling

def evaluate(self, command=None):
    """Evaluate a command and request user input for next command. If no command, return. This function is used for debugging a System object."""    
    # Done evaluating if no command, exit debugger if 'exit'
    if command is None:
        Next = CS.next('Next: ') + f'{repr(self)}\n'
        info = CS.info("Enter to continue or type to evaluate:\n")
        command = input(Next + info + ">>> ")
    
    if command == '':
        return 
    elif command == 'exit':
        raise Stop
    
    if command:
        # Build locals dictionary for evaluating command
        lcs = {} 
        for attr in ('stream', 'unit', 'system'):
            dct = getattr(find, attr)
            lcs.update(dct)
        for key in lcs.keys():
            lcs[key] = lcs[key]()
        lcs['find'] = find
        try:
            out = eval(command, {}, lcs)            
        except Exception as err:
            # Print exception and ask to raise error or continue evaluating
            err = CS.exception(f'{type(err).__name__}:') + f' {str(err)}\n\n'
            info = CS.info(f"Enter to raise error or type to evaluate:\n")
            command = input(err + info + ">>> ")
            if command == '':
                raise err
            evaluate(self, command)        
        else:
            # If successful, continue evaluating
            if out is None:
                if '.show(' not in command:
                    print('\n')
            else:
                print(out)
            command = input(">>> ")
            evaluate(self, command)

def method_debug(self, func):
    """Method decorator for debugging system."""
    def wrapper(*args, **kwargs):
        # Run method and ask to evaluate
        evaluate(self)
        func(*args, **kwargs)
        
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


def notify_run_wrapper(self, func):
    """Decorate a System run method to notify you after each loop"""
    def wrapper(*args, **kwargs):
        if self.recycle:
            func(*args, **kwargs)
            x = self.solver_error['iter']
            input(f'        Finished loop #{x}\n')
        else:
            func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper
    
# %%
    
class _systemUnit(Unit):
    """Dummy unit for displaying a system."""
    line = 'System'
    ID = None
    
_systemUnit._graphics.node['peripheries'] = '2'
_systemUnit._instances = None
del find.line['System']


# %% Process flow

class System:
    """Create a System object that can iteratively run each element in a network of BioSTREAM objects until the recycle stream is converged. A network can have function, Unit and/or System objects. When the network contains an inner System object, it converges/solves it in each loop/iteration.

    **Parameters**

         **ID:** [str] A unique identification

         **network:** tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.

         **recycle:** [Stream] A tear stream for the recycle loop.
         
         **facilities:** tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.

    """
    ### Instance attributes ###
    
    _recycle = None  # Tear stream

    ### Class attributes ###
    
    #: [dict] Dictionary of convergence options regarding maximum number of iterations and molar flow rate and temperature tolerances
    options = SmartBook(units={'Maximum iteration': '#',
                               'Molar tolerance': 'kmol/hr',
                               'Temperature tolerance': 'K'},                  
                        **{'Maximum iteration': 100, 
                           'Molar tolerance': 0.55, 
                           'Temperature tolerance': 0.55})
    
    # [float] Error of the spec objective function
    _spec_error = None

    # [dict] Cached downstream systems by (system, unit, with_facilities) keys
    _cached_downstream_systems = {} 
    
    # [dict] Cached downstream networks by (system, unit) keys
    _cached_downstream_networks = {}

    def __init__(self, ID, network=(), recycle=None, facilities=()):
        #: [dict] Current molar flow and temperature errors and number of iterations made
        self.solver_error = {'mol_error': 0,
                             'T_error': 0,
                             'spec_error': 0,
                             'iter': 0}
        #: set[Unit] All units within the network
        self.units = set()
        
        #: set[Stream] All streams within the network
        self.streams = set()

        #: set[Unit] All units within facilities
        self.offsite_units = set()
        
        #: set[Unit] All streams within facilities
        self.offsite_streams = set()
        
        #: set[System] All subsystems within facilities
        self.offsite_subsystems = set()

        #: set[System] All subsystems in the network
        self.subsystems = set()
        
        self.ID = ID
        self._set_network(network)
        self.recycle = recycle
        self._set_facilities(facilities)
    
    def _set_network(self, network):
        """Set network and cache units, streams, subsystems, feeds and products."""
        # Sort network into surface_units and systems
        units = self.units; units.clear()
        streams = self.streams; streams.clear()
        subsystems = self.subsystems; subsystems.clear()
        for i in network:
            if isinstance(i, Unit):
                units.add(i)
                streams.update(i._ins)
                streams.update(i._outs)
            elif isinstance(i, System):
                subsystems.add(i)
                units.update(i.units)
                streams.update(i.streams)
            elif not isinstance(i, function):
                raise ValueError(f"Only Unit, System, and function objects are valid network elements, not '{type(i).__name__}'")
        streams.discard(missing_stream)        
        
        # Stamp unit in a loop
        if self.recycle:
            for u in units:
                u._in_loop = True
        
        #: set[Stream] All feed streams in the system.
        self.feeds = set(filter(lambda s: not s.source[0] and s.sink[0],
                                   streams))
        
        #: set[Stream] All product streams in the system.
        self.products = set(filter(lambda s: s.source[0] and not s.sink[0],
                                     streams)) 
        
        #: tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.
        self.network = tuple(network)
    
    def _set_facilities(self, facilities):
        """Set facilities and cache offsite units, streams, subsystems, feeds and products."""
        units = self.offsite_units; units.clear()
        streams = self.offsite_streams; streams.clear()
        subsystems = self.offsite_subsystems; subsystems.clear()
        for i in facilities:
            if isinstance(i, Unit):
                units.add(i)
                streams.update(i._ins + i._outs)
            elif isinstance(i, System):
                units.update(i.units)
                streams.update(i.auxiliary_streams)
                subsystems.add(i)
            elif not isinstance(i, function):
                raise ValueError(f"Only Unit, System, and function objects are valid auxiliary elements, not '{type(i).__name__}'")
        
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.
        self.facilities = tuple(facilities)
        
        self.feeds.update(filter(lambda s: not s.source[0] and s.sink[0],
                                 streams))
        self.products.update(filter(lambda s: s.source[0] and not s.sink[0],
                                    streams)) 
    
    @property
    def ID(self):
        """Identification."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        find.system[ID] = self
        self._ID = ID

    @property
    def recycle(self):
        """A tear stream for the recycle loop"""
        return self._recycle

    @recycle.setter
    def recycle(self, stream):
        if isinstance(stream, Stream):
            self._recycle = stream
        elif stream is None:
            self._converge_method = self._run
        else:
            raise ValueError(f"Recycle stream must be a Stream instance or None, not {type(stream).__name__}.")

    @property
    def converge_method(self):
        """Iterative convergence method ('Wegstein', or 'fixed point')."""
        return self._converge_method.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        if self.recycle is None:
            raise ValueError("Cannot set converge method when no recyle is specified")
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' in method:
            self._converge_method = self._Wegstein
        elif 'fixedpoint' in method:
            self._converge_method = self._fixed_point
        else:
            raise ValueError(f"Only 'Wegstein' and 'fixed point' methods are valid, not '{method}'")

    @property
    def _flattened_network(self):
        flattend = []
        for i in self.network:
            if isinstance(i, Unit):
                flattend.append(i)
            elif isinstance(i, System):
                flattend.extend(i._flattened_network)
        return flattend

    
    def _downstream_network(self, unit):
        """Return a list composed of the `unit` and everything downstream."""
        if unit not in self.units:
            return []
        elif self._recycle:
            return self.network
        unit_found = False
        downstream_units = unit._downstream_units
        units = set()
        network = []
        for i in self.network:
            if unit_found:
                if isinstance(i, System):
                    for u in i.units:
                        if u in downstream_units:
                            network.append(i)
                            units.update(i.units)
                            break
                elif i in downstream_units:
                    network.append(i)
                    units.add(i)
                elif (not isinstance(i, Unit)
                    or i.line == 'Balance'):
                    network.append(i)
            else:
                if unit is i:
                    unit_found = True
                    network.append(unit)
                elif isinstance(i, System):
                    if unit in i.units:
                        unit_found = True   
                        network.append(i)
            if not (downstream_units - units): break
        return network

    def _downstream_system(self, unit):
        """Return a system with a network composed of the `unit` and everything downstream (facilities included)."""
        cached = self._cached_downstream_systems
        system = cached.get((self, unit))
        if system: return system
        network = self._downstream_network(unit)
        if not network:
            unit_not_found = True
            for i in self.facilities:
                if unit is i or (isinstance(i, System) and unit in i.units):
                    pos = self.facilities.index(unit)
                    downstream_facilities = self.facilities[pos:]
                    unit_not_found = False
                    break
            if unit_not_found:
                raise ValueError(f'{unit} not found in system')
        else:
            downstream_facilities = self.facilities
        system = System(f'{type(unit).__name__} {unit} and downstream', network,
                        facilities=downstream_facilities)
        cached[unit] = system
        return system
    
    def block(self, element):
        """Create a Block object that can simulate the element and the system downstream.
        
        **Parameter**
        
            **element:** [Unit or Stream] Element in system network or facilities.
        
        """
        return Block(element, self)
    
    def grid(self, *blockfunc_argspace):
        """Create a Grid object that can simulate the argument space for each block function.
        
        **Parameters**
        
            **blockfunc_argspace:** [(function, args)] Iterable of block fuctions and respective arguments.
        
        """
        return Grid(self, *blockfunc_argspace)
    
    def _minimal_diagram(self):
        """Minimally display the network as a box."""
        outs = []
        ins = []
        for s in self.streams:
            source = s.source[0]
            sink = s.sink[0]
            if source in self.units and sink not in self.units:
                outs.append(s)
            elif sink in self.units and source not in self.units:
                ins.append(s)
        subsystem_unit = _systemUnit(self.ID, outs, ins)
        subsystem_unit.line = 'System'
        subsystem_unit.diagram()
        # Reconnect how it was
        for u in self.units:
            u.ins = u._ins
            u.outs = u._outs

    def _surface_diagram(self):
        """Display only surface elements listed in the network."""
        # Get surface items to make nodes and edges
        units = set()
        streams = set()        
        for i in self.network:
            if isinstance(i, Unit):
                units.add(i)
                streams.update(i._ins)
                streams.update(i._outs)
            elif isinstance(i, System):
                outs = []
                ins = []
                for s in i.streams:
                    source = s.source[0]
                    sink = s.sink[0]
                    if source in i.units and sink not in i.units:
                        streams.add(s)
                        outs.append(s)
                    elif sink in i.units and source not in i.units:
                        streams.add(s)
                        ins.append(s)
                subsystem_unit = _systemUnit(i.ID, outs, ins)
                subsystem_unit.line = 'System'
                units.add(subsystem_unit)
        System('', units)._thorough_diagram()
        # Reconnect how it was
        for u in self.units:
            u.ins = u._ins
            u.outs = u._outs
      
    def _thorough_diagram(self):
        """Thoroughly display every unit within the network."""
        # Create a digraph and set direction left to right
        f = Digraph(format='svg')
        f.attr(rankdir='LR')

        # Set up unit nodes
        U = {}  # Contains units by ID
        UD = {}  # Contains full description (ID and line) by ID
        for unit in self.units:
            graphics = unit._graphics
            if unit.ID == '' or not graphics.in_system:
                continue  # Ignore Unit
            
            # Initialize graphics and make Unit node with attributes
            unit._graphics.node_function(unit)
            Type = graphics.name if graphics.name else unit.line
            name = unit.ID + '\n' + Type
            f.attr('node', **unit._graphics.node)
            f.node(name)
            U[unit.ID] = unit
            UD[unit.ID] = name
            
        keys = UD.keys()

        # Set attributes for graph and streams
        f.attr('node', shape='rarrow', fillcolor='#79dae8', orientation='0',
               style='filled', peripheries='1', color='black')
        f.attr('graph', splines='normal', overlap='orthoyx',
               outputorder='edgesfirst', nodesep='0.15', maxiter='10000')
        f.attr('edge', dir='foward')

        for stream in self.streams:
            if stream.ID == '' or stream.ID == 'Missing Stream':
                continue  # Ignore stream

            oU, oi = stream._source
            dU, di = stream._sink
            if oU: oU = oU.ID
            if dU: dU = dU.ID

            # Make stream nodes / unit-stream edges / unit-unit edges
            if oU not in keys and dU not in keys:
                # Stream is not attached to anything
                continue
            elif oU not in keys:
                # Feed stream case
                f.attr('node', shape='rarrow', fillcolor='#79dae8',
                       style='filled', orientation='0', width='0.6',
                       height='0.6', color='black')
                f.node(stream.ID)
                edge_in = U[dU]._graphics.edge_in
                f.attr('edge', arrowtail='none', arrowhead='none',
                       tailport='e', **edge_in[di])
                f.edge(stream.ID, UD[dU])
            elif dU not in keys:
                # Product stream case
                f.attr('node', shape='rarrow', fillcolor='#79dae8',
                       style='filled', orientation='0', width='0.6',
                       height='0.6', color='black')
                f.node(stream.ID)
                edge_out = U[oU]._graphics.edge_out
                f.attr('edge', arrowtail='none', arrowhead='none',
                       headport='w', **edge_out[oi])
                f.edge(UD[oU], stream.ID)
            else:
                # Process stream case
                edge_in = U[dU]._graphics.edge_in
                edge_out = U[oU]._graphics.edge_out
                f.attr('edge', arrowtail='none', arrowhead='normal',
                       **edge_in[di], **edge_out[oi])
                f.edge(UD[oU], UD[dU], label=stream.ID)

        x = IPython.display.SVG(f.pipe(format='svg'))
        IPython.display.display(x)
        
    def diagram(self, kind='surface'):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the system.
        
        **Parameters**
        
            **kind:** Must be one of the following:
                * **'thorough':** Thoroughly display every unit within the network
                * **'surface':** Display only surface elements listed in the network
                * **'minimal':** Minimally display the network as a box
        
        """
        if kind == 'thorough':
            return self._thorough_diagram()
        elif kind == 'surface':
            return self._surface_diagram()
        elif kind == 'minimal':
            return self._minimal_diagram()
        else:
            raise ValueError(f"kind must be either 'thorough', 'surface', or 'minimal'.")
            

    # Methods for running one iteration of a loop
    def _run(self):
        """Rigorous run each element of the system."""
        for a in self.network:
            if isinstance(a, Unit): a._run()
            elif isinstance(a, System): a._converge()
            elif isinstance(a, function): a()
        self.solver_error['iter'] += 1
    
    # Methods for convering the recycle stream
    def _fixed_point(self):
        """Converge system recycle using inner and outer loops with fixed-point iteration."""
        # Reused attributes
        recycle = self.recycle
        run = self._run
        solver_error = self.solver_error
        maxiter, mol_tol, T_tol = tuple(self.options.values())

        def set_solver_error():
            solver_error['mol_error'] = mol_error
            solver_error['T_error'] = T_error

        # Outer loop
        while True:
            mol_old = copy(recycle.mol)
            T_old = recycle.T
            run()
            mol_error = sum(abs(recycle.mol - mol_old))
            T_error = abs(recycle.T - T_old)
            if T_error < T_tol and mol_error < mol_tol:
                break
            if solver_error['iter'] > maxiter:
                set_solver_error()
                raise SolverError(f'Could not converge'
                                  + self._error_info())

        set_solver_error()

    def _Wegstein(self):
        """Converge the system recycle iteratively using Wegstein's method."""
        # Reused attributes
        recycle = self.recycle
        run = self._run
        maxiter, mol_tol, T_tol = tuple(self.options.values())
        solver_error = self.solver_error

        def set_solver_error():
            solver_error['mol_error'] = mol_error
            solver_error['T_error'] = T_error

        # Prepare variables
        len_ = recycle._Nspecies + 1
        x0 = np.zeros(len_)
        gx0 = np.zeros(len_)
        x1 = np.zeros(len_)
        gx1 = np.zeros(len_)
        ones = np.ones(len_)
        s = np.ones(len_)

        # First run
        x0[:-1] = recycle.mol
        x0[-1] = recycle.T
        run()
        x1[:-1] = gx0[:-1] = recycle.mol
        x1[-1] = gx0[-1] = recycle.T

        # Check convergence
        mol_error = sum(abs(gx0[:-1] - x0[:-1]))
        T_error = abs(gx0[-1] - x0[-1])
        converged = mol_error < mol_tol and T_error < T_tol
        if converged:
            set_solver_error()
            return

        # Outer loop
        while not converged:
            run()
            gx1[:-1] = recycle.mol
            gx1[-1] = recycle.T

            # Check if converged
            mol_error = sum(abs(gx1[:-1] - x1[:-1]))
            T_error = abs(gx1[-1] - x1[-1])
            converged = mol_error < mol_tol and T_error < T_tol

            if solver_error['iter'] > maxiter:
                set_solver_error()
                raise SolverError(f'Could not converge'
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
            recycle.mol = x1[:-1]
            recycle.T = x1[-1]
        set_solver_error()
    
    # Default converge method
    _converge_method = _Wegstein
    
    @notify_error
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
        solver_error = self.solver_error

        def error(val):
            setter(val)
            converge()
            solver_error['spec_error'] = e = getter()
            return e

        def converge_spec():
            converge()
            x = solver(error, **kwargs)
            if solver is newton: kwargs['x0'] = x

        self._converge_method = converge_spec

    def _reset_iter(self):
        self.solver_error['iter'] = 0
        for system in self.subsystems:
            system._reset_iter()
    
    def reset_names(self, unit_format=None, stream_format=None):
        """Reset names of all streams and units according to the network order."""
        Unit._default_ID = unit_format if unit_format else ['U', 1]
        Stream._default_ID = stream_format if stream_format else ['S', 1]
        subsystems = set()
        streams = set()
        units = set()
        for i in self.network:
            if isinstance(i, Unit) and i not in units:
                i.ID = ''
                for s in (i._ins + i._outs):
                    streams.add(s)
                    if (s is not missing_stream
                        and s.sink[0] and s.source[0]
                        and s not in streams):
                        s.ID = ''  
            elif isinstance(i, System) and i not in subsystems:
                subsystems.add(i)
                i.reset_names(Unit._default_ID, Stream._default_ID) 
    
    def reset_flows(self):
        """Reset all process streams to zero flow."""
        self.solver_error = {'mol_error': 0,
                             'T_error': 0,
                             'spec_error': 0,
                             'iter': 0}
        for stream in self.streams:
            if stream.source[0]: stream.empty()

    def simulate(self):
        """Converge the network and simulate all units."""
        self._reset_iter()
        for u in self.units:
            u._setup_linked_streams()
            if u._kwargs != u.kwargs:
                u._setup()
                u._kwargs = copy(u.kwargs)
        self._converge()
        for u in self.units:
            u._summary()
        for i in self.facilities:
            if isinstance(i, function):
                i()
            else:
                i.simulate()
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        self._run = notify_run_wrapper(self, self._run)
        self.network = network = list(self.network)
        for i, item in enumerate(network):
            if isinstance(item, Unit):
                item._run = method_debug(item, item._run)
            elif isinstance(item, System):
                item._converge = method_debug(item, item._converge)
            elif isinstance(item, function):
                network[i] = method_debug(item, item)

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
        return self.ID
    
    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + self.ID + '>'

    def show(self):
        """Print all specifications."""
        print(self._info())

    def _error_info(self):
        """Return information on convergence."""
        x = self.solver_error
        recycle = self.recycle

        error_info = ''

        spec = x['spec_error']
        if spec:
            error_info += f'\n specification error: {spec:.3g}'

        if recycle:
            error_info += f"\n convergence error: Flow rate   {x['mol_error']:.2e} {CS.dim(mol_units)}"
            error_info += f"\n                    Temperature {x['T_error']:.2e} {CS.dim(T_units)}"

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
        i = -1; last_i = 0
        while True:
            i = network.find(', ', i+1)
            if (i-last_i) > 30:
                network = (network[:i] + '%' + network[i:])
                last_i = i
            elif i == -1: break
            i += 5
        network = network.replace('%, ', ',\n'+' '*11)
            
        facilities = strtuple(self.facilities)
        i = -1; last_i = 0
        while True:
            i = facilities.find(', ', i+1)
            if (i - last_i) > 30:
                facilities = (facilities[:i] + '%' + facilities[i:])
                last_i = i
            elif i == -1: break
            i += 5
        facilities = facilities.replace('%, ', ',\n'+' '*14)
        
        return (f"System: {self.ID}"
                + str(recycle) + '\n'
                + f" network: {network}\n"
                +(f" facilities: {facilities}" if self.facilities else '')
                + str(error)[1:])
