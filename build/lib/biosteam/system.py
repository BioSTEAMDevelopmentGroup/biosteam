# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
from copy import copy
import IPython
from scipy.optimize import brentq
from graphviz import Digraph
from biosteam.exceptions import Stop, SolverError, notify_error
from biosteam.find import find
from biosteam.stream import Stream, mol_units, T_units
from biosteam.unit import Unit
from biosteam import np
from biosteam.utils import color_scheme, missing_stream
from bookkeep import SmartBook
CS = color_scheme

# Get the function type
function = type(find)

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
        self.network = network
        self.recycle = recycle
        self.facilities = facilities
    
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
    def network(self):
        """A network of BioSTEAM of objects and/or functions that are run element by element until the recycle converges. See example for details."""
        return self._network

    @network.setter
    def network(self, network):
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
        self._network = tuple(network)

    @property
    def facilities(self):
        return self._facilities
    
    @facilities.setter
    def facilities(self, facilities):
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
        self._facilities = tuple(facilities)
        self.feeds.update(filter(lambda s: not s.source[0] and s.sink[0],
                                 streams))
        self.products.update(filter(lambda s: s.source[0] and not s.sink[0],
                                    streams)) 

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
        subsystems = self.subsystems        
        for i in self.network:
            if isinstance(i, Unit):
                units.add(i)
                streams.update(i._ins)
                streams.update(i._outs)
        for sys in subsystems:
            outs = []
            ins = []
            for s in sys.streams:
                source = s.source[0]
                sink = s.sink[0]
                if source in sys.units and sink not in sys.units:
                    streams.add(s)
                    outs.append(s)
                elif sink in sys.units and source not in sys.units:
                    streams.add(s)
                    ins.append(s)
            subsystem_unit = _systemUnit(sys.ID, outs, ins)
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
        
    def diagram(self, kind='thorough'):
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
        for a in self._network:
            if isinstance(a, Unit):
                a._run()
            elif isinstance(a, System):
                a._converge()
            elif isinstance(a, function):
                a()
        self.solver_error['iter'] += 1
    
    def _run_units(self):
        for a in self._network:
            if isinstance(a, Unit):
                a._run()
            elif isinstance(a, System):
                a._run_units()
    
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
        return self._converge_method()

    def set_spec(self, func, var_setter, var_min, var_max, brentq_kwargs={'xtol': 10**-6}):
        """Wrap a bounded solver around the converge method. The solver varies a variable until the desired spec is satisfied.

        **Parameters**

             func: [function] Objective function which should return 0 once the spec is met.

             var_setter: [function] Changes the independent variable.

             var_min: [function] Returns the minimun value of the independent variable.

             var_max: [function] Returns the maximum value of the independent variable.

             brentq_kwargs: [dict] Key word arguments passed to scipy.optimize.brentq solver.

        """
        converge = self._converge_method
        solver_error = self.solver_error

        def error(val):
            var_setter(val)
            converge()
            solver_error['spec_error'] = e = func()
            return e

        def converge_spec():
            converge()
            brentq(error, var_min(), var_max(), **brentq_kwargs)

        self._converge_method = converge_spec

    def _reset_iter(self):
        self.solver_error['iter'] = 0
        for system in self.subsystems:
            system._reset_iter()
    
    def _reset_names(self, unit_format=None, stream_format=None):
        """Reset names of all streams and units in order of network."""
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
                i._reset_names(Unit._default_ID, Stream._default_ID) 
    
    def reset(self):
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
            if u._kwargs != u.kwargs:
                u._setup()
                u._kwargs = copy(u.kwargs)
        self._converge()
        for u in self.units:
            u._operation()
            u._design()
            u._cost()
            u._summary()
        for i in self._facilities:
            if isinstance(i, function):
                i()
            else:
                i.simulate()
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        self._run = notify_run_wrapper(self, self._run)
        self._network = network = list(self._network)
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
        network = self._network
        for i, item in enumerate(network):
            if isinstance(item, Unit):
                item._run = item._run._original
            elif isinstance(item, System):
                item._converge = item._converge._original
            elif isinstance(item, function):
                network[i] = item._original
        self._network = tuple(network)
    
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
            error_info += f"\n convergence error: mol={x['mol_error']:.2e} ({mol_units}), T={x['T_error']:.2e} ({T_units})"

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
        network = str(tuple(str(n) for n in self._network)).replace("'", "")
        return (f"System: {self.ID}"
                + str(recycle) + '\n'
                + f" network: {network}"
                + str(error))
