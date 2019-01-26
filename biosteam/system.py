# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
import time
import copy
import IPython
import sys as SYS
from scipy.optimize import brentq
from graphviz import Digraph
from biosteam.exceptions import Stop, SolverError
from biosteam.lookup import LookUp
from biosteam.stream import Stream, mol_units, T_units
from biosteam.unit import Unit
from biosteam import np
from biosteam.utils import organized_list, get_streams, color_scheme
CS = color_scheme

# Get the function type
def m(): pass
function = type(m)
del m

# %% Debugging and exception handling

def evaluate(command=None):
    """Evaluate a command and request user input for next command. If no command, return. This function is used for debugging a System object."""
    if command is None:
        command = input(CS.info("Enter to continue or type to evaluate:\n") + CS.request(">>> ")) 
        
    # Done evaluating if no command, exit debugger if 'exit'
    if command == '':
        return 
    elif command == 'exit':
        raise Stop

    try:
        # Build locals dictionary for evaluating command
        lcs = {} 
        for attr in ('stream', 'unit', 'system'):
            dct = getattr(LookUp, attr)
            lcs.update(dct)
        for key in lcs.keys():
            lcs[key] = lcs[key]()
        lcs['LookUp'] = LookUp
        
        # Evaluate command and print
        out = eval(command, {}, lcs)
        if out is not None:
            print(out)
    
    except Exception as err:
        # Print exception and ask to raise error or continue evaluating
        print(CS.exception(f'\n{type(err).__name__}:') + f' {str(err)}')
        command = input(CS.info(f"Enter to raise error or type to evaluate:") + CS.request("\n>>> "))
        if command == '':
            raise err
        evaluate(command)
    
    else:
        # If successful, continue evaluating
        evaluate()

def system_debug(self, func):
    """Unit class run method decorator for debugging system."""
    def wrapper(*args):
        # Run method, show info and ask to evaluate
        print()
        func(*args)
        self.show()
        evaluate()
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


def add_func_debugger(self, func):
    """Decorate instance methods to provide a location summary when an error occurs."""
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            out = str(e) + '\n ' + self._info()
            raise type(e)(out).with_traceback(SYS.exc_info()[2])
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def wait(start=''):
    """Make you wait a second and prints dots."""
    print(start, end='')
    time.sleep(0.3)
    for iter in range(3):
        print('.', end='')
        time.sleep(0.2)


def notify_run_wrapper(self, func):
    """Decorate a System run method to notify you after each loop"""
    def wrapper(*args, **kwargs):
        if self.recycle:
            func(*args, **kwargs)
            x = self.solver_error['iter']
            input(CS.note(f'        Finished loop #{x}\n'))
        else:
            func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


def notify_wrapper(self, func):
    """Decorate a convege method of systems within a system for debugging."""
    def wrapper(*args, **kwargs):
        print(f"\nRunning {type(self).__name__} '{self.ID}'", end='')
        wait()
        func(*args, **kwargs)
        print(f" finished!", end='')
        time.sleep(0.4)
        err_info = self._error_info()
        if err_info:
            print(err_info)
        else:
            print()
        evaluate()
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    return wrapper


# %% Metaclass for System format

class metaSystem(type):
    
    @property
    def converge_method(cls):
        """Iterative convergence method ('Wegstein', or 'fixed point')."""
        return cls.converge.__name__[1:]

    @converge_method.setter
    def converge_method(cls, method):
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' in method:
            cls.converge = cls._Wegstein
        elif 'fixedpoint' in method:
            cls.converge = cls._fixed_point
        else:
            raise ValueError(f"Only 'Wegstein' and 'fixed point' methods are valid, not '{method}'")


# %% Process flow

class System(metaclass=metaSystem):
    """Create a System object that can iteratively run each element in a network of BioSTREAM objects until the recycle stream is converged. A network can have function, Unit and/or System objects. When the network contains an inner System object, it converges/solves it in each loop/iteration.

    **Parameters**

         **ID:** [str] A unique identification

         **network:** iterable[Unit, function and/or System] A network that is run element by element until the recycle converges.

         **recycle:** [Stream] A tear stream for the recycle loop.

    **Class Attributes**

         **maxtier** = 100: [int] Maximum number of iterations for converge method.

         **mol_tol** = 0.55: [float] Absolute material tolerance (kmol/hr) of converge method.

         **T_tol** = 0.55: [float] Absolute temperature tolerance (K) of converge method.

    """
    ### Instance attributes ###
    
    _recycle = None  # Tear stream

    ### Class attributes ###
    
    #: [ColorPalette] Manages syntax coloring for the debug method.
    color_scheme = CS
    
    # [iter] Maximum number of iterations for converge method
    maxiter = 100
    
    # [float] Absolute material tolerance (kmol/hr) of converge method
    mol_tol = 0.55
    
    # [float] Absolute temperature tolerance (K) of converge method
    T_tol = 0.55
    
    # [float] Error of the spec objective function
    _spec_error = None

    def __init__(self, ID, network=(), recycle=None):
        #: [dict] Current molar flow and temperature errors and number of iterations made
        self.solver_error = {'mol_error': 0,
                             'T_error': 0,
                             'spec_error': 0,
                             'iter': 0}
        #: list[Unit] All inner and outer units in the network
        self._units = []

        #: list[Stream] All streams in the network
        self._streams = ()

        #: list[Unit] all outer Unit elements of the network
        self._surface_units = []

        #: list[System] all outer System elements of the network
        self._systems = []

        self.ID = ID
        self.network = network
        self.recycle = recycle

    @property
    def streams(self):
        """All streams in the network."""
        return tuple(self._streams)
    
    @property
    def units(self):
        """All inner and outer units in the network."""
        return self._units

    @property
    def ID(self):
        """Identification."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        LookUp.system[ID] = self
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
            self.converge = self._run
        else:
            raise ValueError(f"Recycle stream must be a Stream instance or None, not {type(stream).__name__}.")

    @property
    def network(self):
        """A network of BioSTEAM of objects and/or functions that are run element by element until the recycle converges. See example for details."""
        return self._network

    @network.setter
    def network(self, network):
        if network is None:
            return
        self._network = tuple(network)
        
        # Sort network into surface_units and systems
        _surface_units = self._surface_units
        _systems = self._systems
        for i in network:
            if isinstance(i, Unit):
                _surface_units.append(i)
            elif isinstance(i, System):
                _systems.append(i)
            elif not isinstance(i, function):
                raise ValueError(f"Only Unit, System, and function objects are valid network elements, not '{type(i).__name__}'")
        self._surface_units = list(set(_surface_units))
        self._systems = _systems
        
        # Set all stream and units
        units = copy.copy(self._surface_units)
        streams = get_streams(units)
        for sys in _systems:
            units += sys._units
            streams += sys._streams
        self._units = units = organized_list(units)
        self._streams = organized_list(streams)
        
        # Stamp unit in a loop
        if self.recycle:
            for u in units:
                u._in_loop = True

    @property
    def converge_method(self):
        """Iterative convergence method ('Wegstein', or 'fixed point')."""
        return self.converge.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' in method:
            self.converge = self._Wegstein
        elif 'fixedpoint' in method:
            self.converge = self._fixed_point
        else:
            raise ValueError(f"Only 'Wegstein' and 'fixed point' methods are valid, not '{method}'")

    @property
    def diagram(self, name='diagram'):
        """A Graphviz graph of the network."""
        # Create a digraph and set direction left to right
        f = Digraph(name=name, filename=name, format='svg')
        f.attr(rankdir='LR')

        # Set up unit nodes
        U = {}  # Contains units by ID
        UD = {}  # Contains full description (ID and line) by ID
        for unit in self.units:
            if unit.ID == '':
                continue  # Ignore Unit
            
            if not unit._Graphics.in_system:
                unit.diagram   # Display separately
                continue
            
            # Fill dictionaries
            line = unit.line.replace('_', '\n').replace(' ', '\n')
            name = unit.ID + '\n' + line
            U[unit.ID] = unit
            UD[unit.ID] = name
            
            # Initialize graphics and make Unit node with attributes
            unit._Graphics.node_function(unit)
            f.attr('node', **unit._Graphics.node)
            f.node(name)
            
        keys = UD.keys()
        streams = self.streams

        # Set attributes for graph and streams
        f.attr('node', shape='rarrow', fillcolor='#79dae8', style='filled')
        f.attr('graph', splines='normal', overlap='orthoyx',
               outputorder='edgesfirst', nodesep='0.15', maxiter='10000')
        f.attr('edge', dir='foward')

        for stream in streams:
            if stream.ID == '' or stream.ID == 'Missing Stream':
                continue  # Ignore stream

            oU, oi = stream._source
            dU, di = stream._sink

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
                edge_in = U[dU]._Graphics.edge_in
                f.attr('edge', arrowtail='none', arrowhead='none',
                       tailport='e', **edge_in[di])
                f.edge(stream.ID, UD[dU])
            elif dU not in keys:
                # Product stream case
                f.attr('node', shape='rarrow', fillcolor='#79dae8',
                       style='filled', orientation='0', width='0.6',
                       height='0.6', color='black')
                f.node(stream.ID)
                edge_out = U[oU]._Graphics.edge_out
                f.attr('edge', arrowtail='none', arrowhead='none',
                       headport='w', **edge_out[oi])
                f.edge(UD[oU], stream.ID)
            else:
                # Process stream case
                edge_in = U[dU]._Graphics.edge_in
                edge_out = U[oU]._Graphics.edge_out
                f.attr('edge', arrowtail='none', arrowhead='normal',
                       **edge_in[di], **edge_out[oi])
                f.edge(UD[oU], UD[dU], label=stream.ID)

        x = IPython.display.SVG(f.pipe(format='svg'))
        IPython.display.display(x)

    # Methods for running one iteration of a loop
    def _run(self):
        """Rigorous run each element of the system."""
        for a in self._network:
            if isinstance(a, Unit):
                a.run()
            elif isinstance(a, System):
                a.converge()
            elif isinstance(a, function):
                a()
        self.solver_error['iter'] += 1
    
    # Methods for convering the recycle stream
    def _fixed_point(self):
        """Converge system recycle using inner and outer loops with fixed-point iteration."""
        # Reused attributes
        recycle = self.recycle
        run = self._run
        solver_error = self.solver_error
        T_tol = self.T_tol
        mol_tol = self.mol_tol
        maxiter = self.maxiter

        def set_solver_error():
            solver_error['mol_error'] = mol_error
            solver_error['T_error'] = T_error

        # Outer loop
        while True:
            mol_old = copy.copy(recycle.mol)
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
        mol_tol = self.mol_tol
        T_tol = self.T_tol
        solver_error = self.solver_error
        maxiter = self.maxiter
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
    
    #: Converge the system recycle using an iterative solver (Wegstein by default). 
    converge = _Wegstein

    # Solving specification
    def set_spec(self, func, var_setter, var_min, var_max, brentq_kwargs={'xtol': 10**-6}):
        """Wrap a bounded solver around the converge method. The solver varies a variable until the desired spec is satisfied.

        **Parameters**

             func: [function] Objective function which should return 0 once the spec is met.

             var_setter: [function] Changes the independent variable.

             var_min: [function] Returns the minimun value of the independent variable.

             var_max: [function] Returns the maximum value of the independent variable.

             brentq_kwargs: [dict] Key word arguments passed to scipy.optimize.brentq solver.

        """
        brentq_ = add_func_debugger(self, brentq)
        converge = self.converge
        solver_error = self.solver_error

        def error(val):
            var_setter(val)
            converge()
            solver_error['spec_error'] = func()
            return solver_error['spec_error']

        def converge_spec():
            converge()
            solver_error['spec_error'] = func()
            if abs(solver_error['spec_error']) < brentq_kwargs['xtol']:
                return
            min_ = var_min()
            max_ = var_max()
            brentq_(error, min_, max_, **brentq_kwargs)

        self.converge = converge_spec

    def run_results(self):
        """Run all design basis algorithms (operation, design, and cost methods) for all units."""
        for unit in self.units:
            unit.operation()
            unit.design()
            unit.cost()
        
    # Debugging
    def _debug_on(self):
        """Turn on debug mode."""
        self._run = notify_run_wrapper(self, self._run)
        for unit in self._surface_units:
            unit.run = system_debug(unit, unit.run)
        for sys in self._systems:
            sys.converge = notify_wrapper(sys, sys.converge)

    def _debug_off(self):
        """Turn off debug mode."""
        self._run = self._run._original
        for unit in self._surface_units:
            unit.run = unit.run._original
        for sys in self._systems:
            sys.converge = sys.converge._original
    
    def debug(self):
        """Convege in debug mode. Just try it!"""
        self._debug_on()
        try:
            evaluate()
            self.converge()
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
            unit, i = self.recycle._source
            recycle = f"\n recycle: {unit}-{i}"
        error = self._error_info()
        network = str(tuple(str(n) for n in self._network)).replace("'", "")
        return (f"System: {self.ID}"
                + str(recycle) + '\n'
                + f" network: {network}"
                + str(error))
