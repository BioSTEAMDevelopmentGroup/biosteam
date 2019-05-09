# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:40:28 2018

@author: yoelr
"""

import re
import os
import numpy as np
import pandas as pd
from IPython import display
from graphviz import Digraph
from ._exceptions import _notify_error, DesignWarning
from ._flowsheet import find
from ._graphics import Graphics, default_graphics
from ._stream import Stream
from ._proxy_stream import ProxyStream
from ._heat_utility import HeatUtility
from ._utils import Ins, Outs, missing_stream
from ._power_utility import PowerUtility
from warnings import warn
import biosteam

dir_path = os.path.dirname(os.path.realpath(__file__)) + '\\'


__all__ = ('Unit',)
_do_nothing = lambda self: None

# %% Bounds checking

def _warning(source, msg, category=Warning):
        """Return a Warning object with source description."""
        if isinstance(source, str):
            msg = f'@{source}: ' + msg
        elif source:
            msg = f'@{type(source).__name__} {str(source)}: ' + msg
        return category(msg)

def _warnbounds(key, value, units, bounds, stacklevel, source):
        """A lower level, functional version of "_checkbounds". Issue a warning if value is out of bounds.
        
        **Parameters**
        
            **key:** [str] Name of value
            
            **value:** [number, Quantity, or array]
            
            **units:** [str] Units of measure
            
            **bounds:** [array or Quantity-array] Lower and upper bounds
            
            **stacklevel:** [int] Stacklevel for warning.
        
            **source:** [str or object] Short description or object it describes.
        
        """
        lb, ub = bounds
        # Warn when value is out of bounds
        if not lb<=value and ub>=value:
            units = ' ' + units if units else ''
            try:
                msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb:.4g} to {ub:.4g}{units})."
            except:  # Handle format errors
                msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb} to {ub}{units})."
            warn(_warning(source, msg, DesignWarning), stacklevel=stacklevel)
            
def _lb_warning(key, value, units, lb, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (minimum {lb}{units})."
    
    warn(_warning(source, msg, DesignWarning), stacklevel=stacklevel)
    
def _ub_warning(key, value, units, ub, stacklevel, source):
    units = ' ' + units if units else ''
    try:
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub:.4g}{units})."
    except:  # Handle format errors
        msg = f"{key} ({value:.4g}{units}) is out of bounds (maximum {ub}{units})."
    
    warn(_warning(source, msg, DesignWarning), stacklevel=stacklevel)

def _checkbounds(self, key, value):
        """Issue a warning if value is out of bounds.
        
        **Parameters**
        
            **key:** [str] Name of value
            
            **value:** [number, Quantity, or array]
            
        """
        bounds = self._bounds.get(key)
        if bounds is not None:
            _warnbounds(key, value, self._units.get(key, ''), bounds, 4, self)

def _boundsignore(*args, **kwargs): pass


# %% Unit metaclass

_Unit_is_done = False

class metaUnit(type):
    """Unit metaclass for wrapping up methods with error notifiers, adding key word arguments, and keeping track for Unit lines and inheritance."""
    _enforce_bounds = True
    _CEPCI = 567.5 # Chemical engineering plant cost index (567.5 at 2017)
    def __new__(mcl, clsname, superclasses, new_definitions):
        """Prepare unit methods with wrappers for error notification, and add _kwargs as key word arguments to __init__. """

        if not _Unit_is_done:
            # Abstract Unit class
            cls = type.__new__(mcl, clsname, superclasses, new_definitions)
        else:
            # Add error notification wrapper to every new unit method
            for method_name in ('_setup', '_run', '_operation', '_design', '_cost'):
                if new_definitions.get(method_name):
                    new_definitions[method_name] = _notify_error(new_definitions[method_name])
    
            # Make new Unit class
            cls = type.__new__(mcl, clsname, superclasses, new_definitions)
            
            if cls.__doc__ is Unit.__doc__:
                # Do not inherit docstring from Unit
                cls.__doc__ = None
            if cls.__doc__:
                cls.__doc__ = cls.__doc__.replace('**Parameters**', 
            "**Parameters**\n\n" +
            
        "        **ID:** [str] Unique identification. If set as '', a default ID will be chosen.\n\n" +

        "        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.\n\n" +
        
        "        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.")
    
            # Set line
            # default_line constitutes a new Unit class
            line = cls.line
            if line is 'Unit':
                line = cls.__name__.replace('_', ' ')
                # Set new graphics object for new line
                if not new_definitions.get('_graphics'):
                    cls._graphics = Graphics()
            elif new_definitions.get('line'):
                # Set new graphics for specified line
                if not new_definitions.get('_graphics'):
                    cls._graphics = Graphics()
            
            cls.line = line = re.sub(r"\B([A-Z])", r" \1", line).capitalize()
            kwargs = new_definitions.get('_kwargs')
            if kwargs and '__init__' not in new_definitions:
                # Key word arguments to replace
                inputs = ''
                for key in kwargs.keys():
                    inputs += ', ' + key + '=' + key
        
                # Begin changing __init__ to have kwargs by making a string to execute
                str2exec = f"def __init__(self, ID='', outs=(), ins=None{inputs}):\n"
                str2exec+= f"     _(self, ID, outs, ins{inputs})"
        
                # Execute string and replace __init__
                globs = {'_': Unit.__init__}
                globs.update(kwargs)
                locs = {}
                exec(str2exec, globs, locs)
                cls.__init__ = locs['__init__']
        return cls
    
    @property
    def enforce_bounds(cls):
        """[bool] True if bounds are checked for all instances and False otherwise."""
        return cls._enforce_bounds
    @enforce_bounds.setter
    def enforce_bounds(cls, val):
        val = bool(val)
        cls._checkbounds = _checkbounds if val else _boundsignore
        cls._enforce_bounds = val
    
    @property
    def CEPCI(cls):
        """Chemical engineering plant cost index (CEPCI)"""
        return cls._CEPCI
    @CEPCI.setter
    def CEPCI(cls, CEPCI):
        cls._CEPCI = CEPCI
    
    def __repr__(cls):
        if cls.line == 'Unit' and cls.__name__ == 'Unit':
            return f'biosteam.{cls.__name__}'
        elif cls.line == cls.__name__ or cls.line == 'Unit':
            return f'Unit.{cls.__name__}'
        else:
            return f'{cls.line}.{cls.__name__}'


# %% Unit Operation


class Unit(metaclass=metaUnit):
    """Abstract parent class for Unit objects. Child objects must contain _setup, _run, _design and _cost methods to setup internal objects, estimate stream outputs of a Unit and find design and cost information. These methods should store information in the '_results' dictionary (an instance attribute).  

    **Parameters**

        **ID:** [str] A unique identification. If ID is an empty string (i.e. '' ), a default ID will be chosen. If ID is None, unit will not be registered in flowsheet.

        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.
        
        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.

        ****kwargs:** Keyword arguments that are stored to instance attribute `_kwargs` and later accessed by _setup, _run, _design, and _cost methods.

    **Class Definitions** 
    
        **line** = [Defaults to the class name of the first child class]: [str] Name denoting the type of Unit class    
    
        **CEPCI** = 567.5: [float] Chemical Engineering Plant Cost Index

        **_kwargs** = {}: [dict] Default keyword arguments

        **_init():**
            Initialize components.
        
        **_setup()**
            Create cached data from kwargs. 

        **_run()**
            Run simulation and update output streams.

        **_prepare()**
            Prepare for design (called before `design`).

        **_design()**
            Add design requirements to results "Design" dictionary.

        **_cost()**
            Add itemized purchse costs to results "Cost" dictionary.
            
        **_spec()**
            Calculate specification cost factors and update purchase prices.
            
        **_join():**
            Join additional purchase prices and utility costs.
        
        **_units** = {}: [dict] Default units for results Operation and Design
        
        **_bounds** = {} [dict] Values should be tuples with lower and upper bounds.
        
        **_N_ins** = 1: [int] Expected number of input streams

        **_N_outs** = 2: [int] Expected number of output streams

        **_N_heat_utilities** = 0: [int] Number of heat utilities  

        **_has_power_utility** = False: [bool] If True, a PowerUtility object is created for every instance.
        
        **_has_cost** = True: [bool] Should be True if it has any associated cost
    
        **_linkedstreams = True** = False: [bool] True if outs are proxy streams linked to ins
    
        **_graphics** = default_graphics: [biosteam Graphics] a Graphics object for diagram representation.
        
    .. Note::
        
       Class argument `_init` is called only once when a Unit object is initialized, while `_setup` is run in the `reset` method to update cached data. The `_run` method is called during recycle loop convergence. The rest of the methods are called in the given order for generating results.

    **ins**
        
        list of input streams
        
    **outs**
    
        list of output streams
    
    **Examples**
    
        :doc:`Creating a Unit`
        
        :doc:`Using -pipe- notation`
        
        :doc:`Inheriting from Unit`
    
    """ 
    ### Abstract Attributes ###
    
    # [dict] Default units for results Operation and Design
    _units = {}
    
    # [bool] Should be True if it has any associated cost
    _has_cost = True
    
    # [int] Expected number of input streams
    _N_ins = 1  
    
    # [int] Expected number of output streams
    _N_outs = 2  
    
    # [bool] True if outs are proxy streams linked to ins
    _linkedstreams = False
    
    # [dict] Values should be tuples with lower and upper bounds for results dictionary.
    _bounds = {}
    
    # [int] number of heat utilities
    _N_heat_utilities = 0
    
    # [PowerUtility] A PowerUtility object, if required
    _power_utility = None
    
    # [bool] If True, a PowerUtility object is created for every instance.
    _has_power_utility = False 
    
    # [dict] default key word arguments that are accessed by simulation method.
    _kwargs = {}
    
    # [biosteam Graphics] a Graphics object for diagram representation.
    _graphics = default_graphics
    
    #: [str] The general type of unit, regardless of class
    line = 'Unit'

    ### Other defaults ###

    # [list] Default ID starting letter and number
    _default_ID = ['U', 1]
    
    # Default ID
    _ID = None 
    
    #: [list] HeatUtility objects associated to unit
    _heat_utilities = None
    
    ### Initialize ###
    
    def __init__(self, ID='', outs=(), ins=None, **kwargs):
        self.ID = ID
        self._kwargs = kwargs
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_results()
        self._init_heat_utils()
        self._init_power_util()
        self._init()
        self._setup()
        self._install()

    def reset(self, **kwargs):
        """Reset unit with new key word arguments."""
        self._kwargs.update(kwargs)
        self._setup()

    def _init_ins(self, ins):
        """Initialize input streams."""
        if ins is None:
            self._ins = Ins(self, (missing_stream for i in range(self._N_ins)))
        elif isinstance(ins, Stream):
            self._ins = Ins(self, ins)
        elif isinstance(ins, str):
            self._ins = Ins(self, Stream(ins))
        elif not ins:
            self._ins = Ins(self, (Stream('') for i in range(self._N_ins)))
        else:
            self._ins = Ins(self, (i if isinstance(i, Stream) else Stream(i) for i in ins))
    
    def _init_outs(self, outs):
        """Initialize output streams."""
        if self._linkedstreams:
            if outs is None:
                self._outs = Outs(self, (ProxyStream('*')
                                  for i in self._N_outs))
            elif not outs:
                self._outs = Outs(self, (ProxyStream('') for i in self._ins))
            elif isinstance(outs, Stream):
                self._outs = Outs(self, ProxyStream.asproxy(outs))
            elif isinstance(outs, str):
                self._outs = Outs(self, ProxyStream(outs))
            else:
                self._outs = Outs(self, (ProxyStream.asproxy(o) if isinstance(o, Stream)
                                         else ProxyStream(o)
                                         for o in outs))
        else:
            if outs is None:
                self._outs = Outs(self, (missing_stream for i in range(self._N_outs)))
            elif not outs:
                self._outs = Outs(self, (Stream('') for i in range(self._N_outs)))
            elif isinstance(outs, Stream):
                self._outs = Outs(self, outs)
            elif isinstance(outs, str):
                self._outs = Outs(self, Stream(outs))
            else:
                self._outs = Outs(self, (i if isinstance(i, Stream) else Stream(i) for i in outs))        
    
    def _init_results(self):
        """Initialize results attribute."""
        self._results = {'Design': {},
                         'Cost': {}}
        self._totalcosts = [0, 0] #: [list] Purchase price (USD) and utility cost (USD/hr)
    
    def _init_heat_utils(self):
        """Initialize heat utilities."""
        if self._N_heat_utilities: 
            self._heat_utilities = [HeatUtility() for i in
                                    range(self._N_heat_utilities)]
        
    def _init_power_util(self):
        """Initialize power utility."""
        if self._has_power_utility:
            self._power_utility = PowerUtility()
    
    def _install(self):
        """Cache objects and/or replace methods for computational efficiency."""
        if not self._has_cost:
            self._summary = Unit._cost
            self.simulate = self._run
        else:
            # Result dictionaries for all utilities
            self._utils = utils = []
            heat_utilities = self._heat_utilities
            power_utility = self._power_utility
            if heat_utilities: utils.extend(heat_utilities)
            if power_utility: utils.append(power_utility)
                
            # Itemized purchase costs
            self._purchase_costs = self._results['Cost'].values()
    
    def _link_streams(self):
        """Setup ProxyStream objects if any."""
        if self._linkedstreams: 
            for i, o in zip(self._ins, self._outs): o.link = i
    
    # Forward pipping
    def __sub__(self, other):
        """Source streams."""
        if isinstance(other, Unit):
            other.ins = self.outs
            return other
        elif type(other) is int:
            return self.outs[other]
        elif isinstance(other, Stream):
            self.outs = other
            return self
        elif isinstance(other, (tuple, list, np.ndarray)):
            if isinstance(other[0], int):
                return [self.outs[i] for i in other]
            else:
                self.outs = other
                return self
        else:
            return other.__rsub__(self)

    def __rsub__(self, other):
        """Sink streams."""
        if type(other) is int:
            return self.ins[other]
        elif isinstance(other, Stream):
            self.ins = other
            return self
        elif isinstance(other, (tuple, list, np.ndarray)):
            if all(isinstance(i, int) for i in other):
                return [self.ins[i] for i in other]
            else:
                self.ins = other
                return self
        else:
            return other.__sub__(self)

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    # Abstract methods
    _init     = _do_nothing
    _setup    = _do_nothing
    _run      = _do_nothing
    _prepare  = _do_nothing
    _design   = _do_nothing
    _finalize = _do_nothing
    _cost     = _do_nothing
    _join     = _do_nothing
    _spec     = _do_nothing
    
    # Summary
    def _summary(self):
        """Calculate all results from unit run."""
        self._prepare()
        self._design()
        self._finalize()

    def _finalize(self):
        """Run all cost methods and finalize purchase and utility cost."""
        self._cost()
        self._spec()
        self._join()
        self._totalcosts[:] = sum(self._purchase_costs), sum(i.cost for i in self._utils)

    def _update_purchase_cost(self):
        self._totalcosts[0] = sum(self._purchase_costs)
    
    def _update_utility_cost(self):
        self._totalcosts[1] = sum(i.cost for i in self._utils)

    def simulate(self):
        """Run rigourous simulation and determine all design requirements."""
        self._link_streams()
        self._run()
        self._summary()

    def results(self, with_units=True):
        """Return key results from simulation as a DataFrame if `with_units` is True or as a Series otherwise."""
        results = self._results
        ID = self.ID
        keys = []; addkey = keys.append
        vals = []; addval = vals.append
        include_GHG = hasattr(self, '_totalGHG')
        if with_units:
            if self._power_utility:
                i = self._power_utility
                addkey(('Power', 'Rate'))
                addkey(('Power', 'Cost'))
                addval(('kW', i.rate))
                addval(('USD/hr', i.cost))
            if self._heat_utilities:
                for i in self._heat_utilities:
                    addkey((i.ID, 'Duty'))
                    addkey((i.ID, 'Flow'))
                    addkey((i.ID, 'Cost'))
                    addval(('kJ/hr', i.duty))
                    addval(('kg/hr', i.flow))
                    addval(('USD/hr', i.cost))
            units = self._units
            results = self._results.copy()
            Cost = results.pop('Cost')
            if include_GHG:
                GHG = results.pop('GHG') if 'GHG' in results else None
            for ko, vo in results.items():
                for ki, vi in vo.items():
                    addkey((ko, ki))
                    addval((units.get(ki, ''), vi))
            for ki, vi in Cost.items():
                addkey(('Cost', ki))
                addval(('USD', vi))
            capital, utility = self._totalcosts
            addkey(('Purchase cost', ''))
            addval(('USD', capital))
            addkey(('Utility cost', ''))
            addval(('USD/hr', utility))
            if include_GHG:
                a, b = self._totalGHG
                GHG_units =  self._GHG_units
                for ko, vo in GHG.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval((GHG_units.get(ko, ''), vi))
                a_key, b_key = GHG_units.keys()
                a_unit, b_unit = GHG_units.values()
                addkey(('Total ' + a_key, ''))
                addval((a_unit, a))
                addkey(('Total ' + b_key, ''))
                addval((b_unit, b))
            
            df = pd.DataFrame(vals,
                              pd.MultiIndex.from_tuples(keys),
                              ('Units', ID))
            df.columns.name = self.line
            return df
        else:
            if self._power_utility:
                i = self._power_utility
                addkey(('Power', 'Rate'))
                addkey(('Power', 'Cost'))
                addval(i.rate)
                addval(i.cost)
            if self._heat_utilities:
                for i in self._heat_utilities:
                    addkey((i.ID, 'Duty'))
                    addkey((i.ID, 'Flow'))
                    addkey((i.ID, 'Cost'))
                    addval(i.duty)
                    addval(i.flow)
                    addval(i.cost)
            if include_GHG:
                GHG = results.pop('GHG') if 'GHG' in results else None
            for ko, vo in results.items():
                for ki, vi in vo.items():
                    addkey((ko, ki))
                    addval(vi)
            capital, utility = self._totalcosts
            addkey(('Purchase cost', ''))
            addval(capital)
            addkey(('Utility cost', ''))
            addval(utility)
            if include_GHG:
                GHG_units =  self._GHG_units
                for ko, vo in GHG.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval(vi)
                a, b = self._totalGHG
                a_key, b_key = GHG_units.keys()
                addkey(('Total ' + a_key, ''))
                addval(a)
                addkey(('Total ' + b_key, ''))
                addval(b)
            series = pd.Series(vals, pd.MultiIndex.from_tuples(keys))
            series.name = ID
            return series

    _checkbounds = _checkbounds
    def _lb_warning(self, key, value, lb):
        """Warn that value is below lower bound."""
        _lb_warning(key, value, self._units.get(key), lb, 4, self)

    def _lb_check(self, key, value):
        """Warn that value if below lower bound."""
        lb = self._bounds[key]
        if value < lb: _lb_warning(key, value, self._units.get(key), 4, self)

    @property
    def CEPCI(self):
        """Chemical engineering plant cost index (CEPCI)."""
        return type(self).CEPCI
    @CEPCI.setter
    def CEPCI(self, CEPCI):
        raise AttributeError('Cannot change class attribute through an instance.')

    @property
    def ID(self):
        """Unique Identification (str). If set as '', it will choose a default ID."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        unit = find.unit
        ID_old = self._ID
        if ID_old and ID_old in unit: del unit[ID_old]
        if ID == '':
            # Select a default ID if requested
            letter, number = self._default_ID
            self._default_ID[1] += 1
            ID = letter + str(number)
            self._ID = ID
            unit[ID] = self
            return
        elif ID:
            ID = ID.replace(' ', '_')
            ID_words = ID.split('_')
            if not all(word.isalnum() for word in ID_words):
                raise ValueError('ID cannot have any special characters')
            unit[ID] = self
            self._ID = ID

    # Input and output streams
    @property
    def ins(self):
        # list of input streams
        return self._ins
    @ins.setter
    def ins(self, streams):
        self._ins.__init__(self, streams)
        
    @property
    def outs(self):
        # list of output streams
        return self._outs
    @outs.setter
    def outs(self, streams):
        self._outs.__init__(self, streams)

    @property
    def _upstream_neighbors(self):
        """Return set of upsteam neighboring units."""
        neighbors = set()
        for s in self._ins:
            u_source = s.source
            if u_source: neighbors.add(u_source)
        return neighbors

    @property
    def _downstream_neighbors(self):
        """Return set of downstream neighboring units."""
        neighbors = set()
        for s in self._outs:
            u_sink = s.sink
            if u_sink: neighbors.add(u_sink)
        return neighbors

    @property
    def _downstream_units(self):
        """Return all units downstreasm."""
        downstream_units = set()
        outer_periphery = self._downstream_neighbors
        inner_periphery = None
        old_length = -1
        new_length = 0
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            downstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                outer_periphery.update(unit._downstream_neighbors)
            new_length = len(downstream_units)
            
        return downstream_units

    def _neighborhood(self, radius=1):
        """Return all neighboring units within given radius.
        
        **Parameters**
        
            **radius:**[int] Maxium number streams between neighbors.
        
        """
        radius -= 1
        if radius < 0: return set()
        upstream_neighbors = self._upstream_neighbors
        downstream_neighbors = self._downstream_neighbors
        neighborhood = upstream_neighbors.union(downstream_neighbors)
        direct_neighborhood = neighborhood
        
        for i in range(radius):
            neighbors = set()
            for neighbor in direct_neighborhood:
                neighbors.update(neighbor._upstream_neighbors)
                neighbors.update(neighbor._downstream_neighbors)
            if neighbors == direct_neighborhood: break
            direct_neighborhood = neighbors
            neighborhood.update(direct_neighborhood)
        
        return neighborhood

    def diagram(self, radius=0, file=None):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the unit and all neighboring units within given radius.
        
        **Parameters**
        
            **radius:** [int] Maxium number streams between neighbors.
        
            **file:** Must be one of the following:
                * [str] File name to save diagram. If format not included, saves file as svg.
                * [None] Display diagram in console.
        
        """
        if radius > 0:
            neighborhood = self._neighborhood(radius)
            neighborhood.add(self)
            sys = biosteam.System('', neighborhood)
            return sys.diagram()
        
        graphics = self._graphics

        # Make a Digraph handle
        f = Digraph(name='unit', filename='unit', format='svg')
        f.attr('graph', ratio='0.5', splines='normal', outputorder='edgesfirst',
               nodesep='1.1', ranksep='0.8', maxiter='1000')  # Specifications
        f.attr(rankdir='LR')  # Left to right

        # If many streams, keep streams close
        if (len(self.ins) >= 3) or (len(self.outs) >= 3):
            f.attr('graph', nodesep='0.4')

        # Initialize node arguments based on unit and make node
        type_ = graphics.node_function(self) or self.line
        name = self.ID + '\n' + type_
        f.attr('node', **self._graphics.node)
        f.node(name)

        # Set stream node attributes
        f.attr('node', shape='rarrow', fillcolor='#79dae8',
               style='filled', orientation='0', width='0.6',
               height='0.6', color='black', peripheries='1')

        # Make nodes and edges for input streams
        di = 0  # Destination position of stream
        for stream in self.ins:
            if not stream.ID or stream.ID == 'Missing Stream': continue
            f.node(stream.ID)  
            edge_in = self._graphics.edge_in
            if di >= len(edge_in): di = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])
            f.edge(stream.ID, name)
            di += 1

        # Make nodes and edges for output streams
        oi = 0  # Origin position of stream
        for stream in self.outs:
            f.node(stream.ID) 
            edge_out = self._graphics.edge_out  
            if oi >= len(edge_out): oi = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(name, stream.ID)
            oi += 1

        # Display digraph on console
        display.display(display.Image(f.pipe(format='png')))
    
    ### Net input and output flows ###
    
    # Molar flow rates
    @property
    def _mol_in(self):
        """Molar flows going in (kmol/hr)."""
        return sum(s.mol for s in self.ins)

    @property
    def _mol_out(self):
        """Molar flows going out (kmol/hr)."""
        return sum(s.mol for s in self.outs)

    @property
    def _molfrac_in(self):
        """Molar fractions going in (kmol/hr)."""
        return self._mol_in/self._molnet_in

    @property
    def _molfrac_out(self):
        """Molar fractions going in."""
        return self._mol_out/self._molnet_out

    @property
    def _molnet_in(self):
        """Net molar flow going in (kmol/hr)."""
        return sum(s.molnet for s in self.ins)

    @property
    def _molnet_out(self):
        """Net molar flow going out (kmol/hr)."""
        return sum(s.molnet for s in self.outs)

    # Mass flow rates
    @property
    def _mass_in(self):
        """Mass flows going in (kg/hr)."""
        return sum(s.mass for s in self.ins)

    @property
    def _mass_out(self):
        """Mass flows going out (kg/hr)."""
        return sum(s.mass for s in self.outs)

    @property
    def _massfrac_in(self):
        """Mass fractions going in."""
        return self._mass_in/self._massnet_in

    @property
    def _massfrac_out(self):
        """Mass fractions going out."""
        return self._mass_out/self._massnet_out

    @property
    def _massnet_in(self):
        """Net mass flow going in (kg/hr)."""
        return sum(s.massnet for s in self.ins)

    @property
    def _massnet_out(self):
        """Net mass flow going out (kg/hr)."""
        return sum(s.massnet for s in self.outs)

    # Volumetric flow rates
    @property
    def _vol_in(self):
        """Volumetric flows going in (m3/hr)."""
        return sum(s.vol for s in self.ins)

    @property
    def _volnet_in(self):
        """Net volumetric flow going in (m3/hr)."""
        return sum(self._vol_in)

    @property
    def _volfrac_in(self):
        """Volumetric fractions going in."""
        return self._vol_in/self._volnet_in

    @property
    def _vol_out(self):
        """Volumetric flows going out (m3/hr)."""
        return sum(s.vol for s in self.outs)

    @property
    def _volnet_out(self):
        """Net volumetric flow going out (m3/hr)."""
        return sum(self._vol_out)

    @property
    def _volfrac_out(self):
        """Volumetric fractions going out."""
        return self._vol_out/self._volnet_out

    # Enthalpy flow rates
    @property
    def _H_in(self):
        """Enthalpy flow going in (kJ/hr)."""
        return sum(s.H for s in self.ins)

    @property
    def _H_out(self):
        """Enthalpy flow going out (kJ/hr)."""
        return sum(s.H for s in self.outs)

    @property
    def _Hf_in(self):
        """Enthalpy of formation flow going in (kJ/hr)."""
        _Hf_in = 0
        for stream in self.ins:
            _Hf_in += stream.Hf
        return _Hf_in

    @property
    def _Hf_out(self):
        """Enthalpy of formation flow going out (kJ/hr)."""
        _Hf_out = 0
        for stream in self.outs:
            _Hf_out += stream.Hf
        return _Hf_out

    @property
    def _Hnet(self):
        """Net enthalpy flow (including enthalpies of formation)."""
        return self._H_out - self._H_in + self._Hf_out - self._Hf_in
    
    # Representation
    def _info(self, **show_units):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        info+= f'ins...\n'
        i = 0
        for stream in self.ins:
            stream_info = stream._info(**show_units)
            if 'Missing Stream' in stream_info:
                info += f'[{i}] {stream.ID}\n'
                i += 1
                continue
            unit = stream._source
            index = stream_info.index('\n')
            source_info = f'  from  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + source_info + stream_info[index+1:] + '\n'
            i += 1
        info += f'outs...\n'
        i = 0
        for stream in self._outs:
            stream_info = stream._info(**show_units)
            if 'Missing Stream' in stream_info:
                info += f'[{i}] {stream.ID}\n'
                i += 1
                continue
            unit = stream._sink
            index = stream_info.index('\n')
            sink_info = f'  to  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + sink_info + stream_info[index+1:] + '\n'
            i += 1
        info = info.replace('\n ', '\n    ')
        return info[:-1]

    def show(self, **show_units):
        """Prints information on unit."""
        print(self._info(**show_units))
    
    def __str__(self):
        if self.ID:
            return self.ID
        else:
            return type(self).__name__

    def __repr__(self):
        if self.ID:
            return f'<{type(self).__name__}: {self.ID}>'
        else:
            return f'<{type(self).__name__}>'

_Unit_is_done = True

