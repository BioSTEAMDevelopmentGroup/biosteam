# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 15:04:55 2018

@author: yoelr
"""
from scipy.optimize import newton
from ._exceptions import SolverError, _try_method
from ._flowsheet import find, make_digraph, save_digraph
from ._stream import Stream
from ._facility import Facility
from ._unit import Unit
from ._report import save_report
from ._utils import colors, MissingStream, strtuple, \
                    conditional_wegstein, conditional_aitken
import biosteam as bst

__all__ = ('System',)

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
        lcs = {} 
        for attr in ('stream', 'unit', 'system'):
            dct = getattr(find, attr).__dict__
            lcs.update({i:j() for i, j in dct.items()})
        lcs.update({attr:getattr(bst, attr) for attr in bst.__all__})
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
    
# %% System node for diagram

class _systemUnit(Unit, isabstract=True):
    """Dummy unit for displaying a system."""
    line = 'System'
    ID = None
    _N_ins = _N_outs = 1
    _ID = property(lambda self: self.ID)

_sysgraphics = _systemUnit._graphics
_sysgraphics.edge_in = _sysgraphics.edge_in * 10
_sysgraphics.edge_out = _sysgraphics.edge_out * 15
_sysgraphics.node['peripheries'] = '2'

class _streamUnit(Unit, isabstract=True):
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


class System(metaclass=system):
    """Create a System object that can iteratively run each element in a network of BioSTREAM objects until the recycle stream is converged. A network can have function, Unit and/or System objects. When the network contains an inner System object, it converges/solves it in each loop/iteration.

    **Parameters**

         **ID:** [str] A unique identification. If ID is None, instance will not be registered in flowsheet.

         **network:** tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.

         **recycle:** [Stream] A tear stream for the recycle loop.
         
         **facilities:** tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.

    """
    ### Class attributes ###
    
    #: Maximum number of iterations
    maxiter = 100

    #: Molar tolerance (kmol/hr)
    molar_tolerance = 0.01
    
    #: Temperature tolerance (K)
    T_tolerance = 0.10

    # [dict] Cached downstream systems by (system, unit, with_facilities) keys
    _cached_downstream_systems = {} 

    def __init__(self, ID, network, recycle=None, facilities=()):
        
        #: Molar flow rate error (kmol/hr)
        self._mol_error = 0
        
        #: Temperature error (K)
        self._T_error = 0
        
        #: Specification error
        self._spec_error = 0
        
        #: Number of iterations
        self._iter = 0
         
        #: set[Stream] All streams within the system
        self.streams = streams = set()

        #: set[System] All subsystems in the system
        self.subsystems = subsystems = set()
        
        #: list[Unit] Network of only unit operations
        self._unitnetwork = units = []
        
        #: tuple[Unit, function and/or System] A network that is run element by element until the recycle converges.
        self.network = tuple(network)
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
        streams.discard(MissingStream) 
        
        # link all unit operations with linked streams
        for u in units:
            if inst(u, bst.Static): u.link_streams()
        
        #: set[Unit] All units within the system
        self.units = units = set(units)
        
        #: set[Unit] All units in the network that have costs
        self._network_costunits = costunits = {i for i in units if i._cost}
        
        #: set[Unit] All units that have costs.
        self._costunits = costunits = costunits.copy()
        
        #: tuple[Unit, function, and/or System] Offsite facilities that are simulated only after completing the network simulation.
        self.facilities = tuple(facilities)
        for i in facilities:
            if inst(i, Unit):
                units.add(i)
                streams.update(i._ins + i._outs)
                if i._cost: costunits.add(i)
                if inst(i, Facility): i._system = self
            elif inst(i, System):
                units.update(i.units)
                streams.update(i.streams)
                subsystems.add(i)
                costunits.update(i._costunits)
        
        #: set[Stream] All feed streams in the system.
        self.feeds = set(filter(_isfeed, streams))
        
        #: set[Stream] All product streams in the system.
        self.products = set(filter(_isproduct, streams)) 
        
        #: [TEA] System object for Techno-Economic Analysis.
        self._TEA = None
        
        if recycle is None:
            self._converge = self._run
        else:
            assert isinstance(recycle, Stream), (
             "recycle must be a Stream instance or None, not "
            f"{type(recycle).__name__}")
        self._recycle = recycle
        
        if ID:
            ID = ID.replace(' ', '_')
            ID_words = ID.split('_')
            assert all(word.isalnum() for word in ID_words), ('ID cannot have any'
                                                              'special characters')
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

    @property
    def converge_method(self):
        """Iterative convergence method ('wegstein', 'aitken', or 'fixed point')."""
        return self._converge.__name__[1:]

    @converge_method.setter
    def converge_method(self, method):
        if self.recycle is None:
            raise ValueError("cannot set converge method when no recyle is specified")
        method = method.lower().replace('-', '').replace(' ', '')
        if 'wegstein' == method:
            self._converge = self._wegstein
        elif 'fixedpoint' == method:
            self._converge = self._fixed_point
        elif 'aitken' == method:
            self._converge = self._aitken
        else:
            raise ValueError(f"only 'wegstein', 'aitken', and 'fixed point' methods are valid, not '{method}'")

    
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
            unit_found = False
            inst = isinstance
            for pos, i in enumerate(self.facilities):
                if unit is i or (inst(i, System) and unit in i.units):
                    downstream_facilities = self.facilities[pos:]
                    unit_found = True
                    break
            assert unit_found, f'{unit} not found in system'
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
                    None, feed)
        _streamUnit('\n'.join([i.ID for i in outs]),
                    product, None)
        unit = _systemUnit(self.ID, feed, product)
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
                    units.add(_streamUnit('\n'.join([i.ID for i in feeds]), None, feed))
                    ins.append(feed)
                else: ins += feeds
                
                if len(products) > 1:
                    product = Stream(None)
                    product._ID = ''
                    units.add(_streamUnit('\n'.join([i.ID for i in products]),
                                          product, None))
                    outs.append(product)
                else: outs += products
                
                subsystem_unit = _systemUnit(i.ID, ins, outs)
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
        
    def diagram(self, kind='surface', file=None, format='png'):
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
    def _iter_run(self, mol):
        """Run the system at specified recycle molar flow rate.
        
        **Parameters**
        
            **mol:** [array] Recycle molar flow rates.
            
        **Returns**
        
            **rmol:** [array] New recycle molar flow rates.
            
            **unconverged:** [bool] True if recycle has not converged.
            
        """
        recycle = self.recycle
        rmol = recycle.mol
        rmol[:] = mol
        T = recycle.T
        self._run()
        self._mol_error = abs(mol - recycle.mol).sum()
        self._T_error = abs(T - recycle.T)
        self._iter += 1
        if self._mol_error < self.molar_tolerance and self._T_error < self.T_tolerance:
            return rmol.copy(), False
        elif self._iter > self.maxiter:
            raise SolverError(f'{repr(self)} could not converge' + self._error_info())
        else:
            return rmol.copy(), True
        
    def _run(self):
        """Rigorous run each element of the system."""
        inst = isinstance
        _try = _try_method
        for a in self.network:
            if inst(a, Unit): _try(a._run)
            elif inst(a, System): a._converge()
            else: a() # Assume it is a function
    
    # Methods for convering the recycle stream
    def _fixed_point(self):
        """Converge system recycle using inner and outer loops with fixed-point iteration."""
        r = self.recycle
        rmol = r._mol
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
        conditional_wegstein(self._iter_run, self.recycle._mol.copy())
    
    def _aitken(self):
        """Converge the system recycle iteratively using Aitken's method."""
        conditional_aitken(self._iter_run, self.recycle._mol.copy())
    
    # Default converge method
    _converge = _aitken

    def set_spec(self, getter, setter, solver=newton, **kwargs):
        """Wrap a solver around the converge method.

        **Parameters**

             **getter:** [function] Returns objective value.

             **setter:** [function] Changes independent variable.
             
             **solver:** [function] Solves objective function.

             **kwargs:** [dict] Key word arguments passed to solver.

        """
        converge = self._converge

        def error(val):
            setter(val)
            converge()
            self._spec_error = e = getter()
            return e

        def converge_spec():
            x = solver(error, **kwargs)
            if solver is newton: kwargs['x0'] = x

        self._converge = converge_spec

    def _reset_iter(self):
        self._iter = 0
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
                if (s and s._sink and s._source
                    and s not in streams):
                    s.ID = ''
                    streams.add(s)
            units.add(i)
    
    def reset_flows(self):
        """Reset all process streams to zero flow."""
        self._iter = 0
        self._spec_error = 0
        self._mol_error = 0
        self._T_error = 0
        
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
            elif callable(item):
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
            elif callable(item):
                network[i] = item._original
        self.network = tuple(network)
    
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
    
    def _ipython_display_(self):
        try: self.diagram('minimal')
        except: pass
        self.show()

    def _error_info(self):
        """Return information on convergence."""
        recycle = self.recycle
        error_info = ''
        if self._spec_error:
            error_info += f'\n specification error: {self._spec_error:.3g}'
        if recycle:
            error_info +=(f"\n convergence error: Flow rate   {self._mol_error:.2e} kmol/hr"
                          f"\n                    Temperature {self._T_error:.2e} K")
        if self._spec_error or recycle:
            error_info += f"\n iterations: {self._iter}"
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
