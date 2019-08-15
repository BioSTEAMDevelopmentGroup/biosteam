# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:40:28 2018

@author: yoelr
"""
import re
import numpy as np
import pandas as pd
from graphviz import Digraph
from ._exceptions import DesignWarning, _try_method
from ._flowsheet import find, save_digraph
from ._graphics import Graphics, default_graphics
from ._stream import Stream
from ._heat_utility import HeatUtility
from ._utils import Ins, Outs, MissingStream, NotImplementedMethod
from ._power_utility import PowerUtility
from warnings import warn
import biosteam as bst

__all__ = ('Unit',)

# %% Unit neighbors

def _add_upstream_neighbors(unit, set):
    """Add upsteam neighboring units to set."""
    for s in unit._ins:
        u_source = s._source
        if u_source: set.add(u_source)

def _add_downstream_neighbors(unit, set):
    """Add downstream neighboring units to set."""
    for s in unit._outs:
        u_sink = s._sink
        if u_sink: set.add(u_sink)


# %% Bounds checking

def _warning(source, msg, category=Warning):
        """Return a Warning object with source description."""
        if isinstance(source, str):
            msg = f'@{source}: ' + msg
        elif source:
            msg = f'@{type(source).__name__} {str(source)}: ' + msg
        return category(msg)
            
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

# %% Metaclass for Unit operations
    
class unit(type):
    """Unit metaclass for keeping track for Unit lines and graphics."""
    def __new__(mcl, name, bases, dct, isabstract=False, unitdoc=True):
        # Make new Unit class
        cls = type.__new__(mcl, name, bases, dct)
        
        try: Unit
        except NameError: return cls
        
        # Set line
        # default_line constitutes a new Unit class
        line = cls.line
        if line in ('Unit', 'Mixer', 'Static', 'Splitter', 'Solids separator', 'Facility'):
            line = cls.__name__.replace('_', ' ')
            # Set new graphics object for new line
            if '_graphics' not in dct: cls._graphics = Graphics.box(cls._N_ins, cls._N_outs)
        elif 'line' in dct and '_graphics' not in dct:
            # Set new graphics for specified line
            cls._graphics = Graphics.box(cls._N_ins, cls._N_outs)
        
        cls.line = re.sub(r"\B([A-Z])", r" \1", line).capitalize()
        
        if unitdoc and cls.__doc__:
            cls.__doc__ = cls.__doc__.replace('**Parameters**', 
   "**Parameters**"
"\n"    
"\n        **ID:** [str] Unique identification. If set as '', a default ID will be chosen."
"\n"
"\n        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given."
"\n"
"\n        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.")
        
        if isabstract: return cls
        
        if not hasattr(cls, '_run'):
            raise NotImplementedError("'Unit' subclass must have a '_run' method unless the 'isabstract' keyword argument is True")
        
        return cls


# %% Unit Operation

class Unit(metaclass=unit):
    """Abstract parent class for Unit objects. Child objects must contain _run, _design and _cost methods to estimate stream outputs of a Unit and find design and cost information.  

    **Parameters**

        **ID:** [str] A unique identification. If ID is an empty string (i.e. '' ), a default ID will be chosen. If ID is None, unit will not be registered in flowsheet.

        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.

        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.
                

    **Class Definitions** 
    
        **line** = [Defaults to the class name of the first child class]: [str] Name denoting the type of Unit class    
        
        **BM** = None: [float] Bare module factor (installation factor)
        
        **_run()**
            Run simulation and update output streams.

        **_design()**
            Add design requirements to "_Design" dictionary attribute.

        **_cost()**
            Add itemized purchse costs to results "_Cost" dictionary attribute.
            
        .. Note::
           
           The `_run` method is called during recycle loop convergence. The `_design` and `_cost` methods are called in the given order for generating results.
        
        **_units** = {}: [dict] Default units for results Operation and Design
        
        **_N_ins** = 1: [int] Expected number of input streams

        **_N_outs** = 2: [int] Expected number of output streams

        **_N_heat_utilities** = 0: [int] Number of heat utilities  

        **_has_power_utility** = False: [bool] If True, a PowerUtility object is created for every instance.
        
        **_has_cost** = True: [bool] Should be True if it has any associated cost
    
        **_graphics** = <Graphics>: [biosteam Graphics] Settings for diagram representation.

    **ins**
        
        list of input streams
        
    **outs**
    
        list of output streams
    
    **Examples**
    
        :doc:`Creating a Unit`
        
        :doc:`Using -pipe- notation`
        
        :doc:`Inheriting from Unit`
        
        :doc:`Unit decorators`
    
    """ 
    ### Abstract Attributes ###
    
    #: [float] Bare module factor (installation factor).
    BM = None
    
    # [dict] Default units for construction
    _units = {}
    
    # [int] Expected number of input streams
    _N_ins = 1  
    
    # [int] Expected number of output streams
    _N_outs = 2  
    
    # [int] number of heat utilities
    _N_heat_utilities = 0
    
    # [PowerUtility] A PowerUtility object, if required
    _power_utility = None
    
    # [bool] If True, a PowerUtility object is created for every instance.
    _has_power_utility = False 
    
    # [biosteam Graphics] a Graphics object for diagram representation.
    _graphics = default_graphics
    
    # [str] The general type of unit, regardless of class
    line = 'Unit'

    ### Other defaults ###

    # [list] Default ID starting letter and number
    _default_ID = ['U', 1]
    
    # Default ID
    _ID = None 
    
    #: [list] HeatUtility objects associated to unit
    _heat_utilities = ()
        
    ### Initialize ###
    
    def __init__(self, ID='', ins=None, outs=()):
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_results()
        self._init_heat_utils()
        self._init_power_util()
        self.ID = ID

    def _init_ins(self, ins):
        """Initialize input streams."""
        if ins is None:
            self._ins = Ins(self, (MissingStream for i in range(self._N_ins)))
        elif isinstance(ins, Stream):
            self._ins = Ins(self, (ins,))
        elif isinstance(ins, str):
            self._ins = Ins(self, (Stream(ins),))
        elif not ins:
            self._ins = Ins(self, (Stream('') for i in range(self._N_ins)))
        else:
            self._ins = Ins(self, (i if isinstance(i, Stream) else Stream(i) for i in ins))
    
    def _init_outs(self, outs):
        """Initialize output streams."""
        if outs is None:
            self._outs = Outs(self, (MissingStream for i in range(self._N_outs)))
        elif not outs:
            self._outs = Outs(self, (Stream('') for i in range(self._N_outs)))
        elif isinstance(outs, Stream):
            self._outs = Outs(self, (outs,))
        elif isinstance(outs, str):
            self._outs = Outs(self, (Stream(outs),))
        else:
            self._outs = Outs(self, (i if isinstance(i, Stream) else Stream(i) for i in outs))        
    
    def _init_results(self):
        """Initialize attributes to store results."""
        # [dict] Updated in `_cost` method
        self._Cost = {}
        
        # [dict] Updated in `_design` method
        self._Design = {}
        
        # [dict] Greenhouse gas emissions
        self._GHGs = {}
    
    def _init_heat_utils(self):
        """Initialize heat utilities."""
        if self._N_heat_utilities: 
            self._heat_utilities = [HeatUtility() for i in
                                    range(self._N_heat_utilities)]
        
    def _init_power_util(self):
        """Initialize power utility."""
        if self._has_power_utility:
            self._power_utility = PowerUtility()         
    
    # Forward pipping
    def __sub__(self, other):
        """Source streams."""
        if isinstance(other, Unit):
            other._ins[:] = self._outs
            return other
        elif type(other) is int:
            return self.outs[other]
        elif isinstance(other, Stream):
            self._outs[:] = (other,)
            return self
        elif isinstance(other, (tuple, list, np.ndarray)):
            if isinstance(other[0], int):
                return [self.outs[i] for i in other]
            else:
                self._outs[:] = other
                return self
        else:
            return other.__rsub__(self)

    def __rsub__(self, other):
        """Sink streams."""
        if type(other) is int:
            return self._ins[other]
        elif isinstance(other, Stream):
            self._ins[:] = (other,)
            return self
        elif isinstance(other, (tuple, list, np.ndarray)):
            if all(isinstance(i, int) for i in other):
                return [self._ins[i] for i in other]
            else:
                self._ins[:] = other
                return self

    # Backwards pipping
    __pow__ = __sub__
    __rpow__ = __rsub__
    
    # Abstract methods
    _design   = NotImplementedMethod
    _cost     = NotImplementedMethod
    
    # Summary
    def _summary(self):
        """Calculate all results from unit run."""
        self._design()
        self._cost()
    
    @property
    def purchase_cost(self):
        """Total purchase cost (USD)."""
        return sum(self._Cost.values())
    
    @property
    def installation_cost(self):
        """Installation cost (USD)."""
        return self.BM * sum(self._Cost.values())
    
    @property
    def utility_cost(self):
        """Total utility cost (USD/hr)."""
        if self._power_utility:
            return (sum([i.cost for i in self._heat_utilities])
                                 + self._power_utility.cost)
        else:
            return sum([i.cost for i in self._heat_utilities])

    def simulate(self):
        """Run rigourous simulation and determine all design requirements."""
        _try_method(self._run)
        _try_method(self._summary)

    def results(self, with_units=True):
        """Return key results from simulation as a DataFrame if `with_units` is True or as a Series otherwise."""
        ID = self.ID
        keys = []; addkey = keys.append
        vals = []; addval = vals.append
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
                    addval(('kmol/hr', i.flow))
                    addval(('USD/hr', i.cost))
            units = self._units
            Cost = self._Cost
            for ki, vi in self._Design.items():
                addkey(('Design', ki))
                addval((units.get(ki, ''), vi))
            for ki, vi in Cost.items():
                addkey(('Cost', ki))
                addval(('USD', vi))
            addkey(('Purchase cost', ''))
            addval(('USD', self.purchase_cost))
            addkey(('Utility cost', ''))
            addval(('USD/hr', self.utility_cost))
            if self._GHGs:
                a, b = self._totalGHG
                GHG_units =  self._GHG_units
                for ko, vo in self._GHGs.items():
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
            for ki, vi in self._Design.items():
                addkey(('Design', ki))
                addval(vi)
            for ki, vi in self._Cost.items():
                addkey(('Cost', ki))
                addval(vi)    
            if self._GHGs:
                GHG_units =  self._GHG_units
                for ko, vo in self._GHGs.items():
                    for ki, vi in vo.items():
                        addkey((ko, ki))
                        addval(vi)
                a, b = self._totalGHG
                a_key, b_key = GHG_units.keys()
                addkey(('Total ' + a_key, ''))
                addval(a)
                addkey(('Total ' + b_key, ''))
                addval(b)
            addkey(('Purchase cost', ''))
            addval(self.purchase_cost)
            addkey(('Utility cost', ''))
            addval(self.utility_cost)
            series = pd.Series(vals, pd.MultiIndex.from_tuples(keys))
            series.name = ID
            return series

    # def __getattr__(self, name):
    #     if name is '_results':
    #         warn(DeprecationWarning("'_results' dictionary will be deprecated in favor of attributes '_Design', '_Cost', and '_GHGs'"))
    #         self._results = results = {'Design': self._Design,
    #                                    'Cost': self._Cost,
    #                                    'GHG': self._GHGs}
    #         return results
    #     elif name in ('_Cost', '_Design') and hasattr(self, '_results'):
    #         warn(DeprecationWarning("'_results' dictionary will be deprecated in favor of attributes '_Design', '_Cost', and '_GHGs'"))
    #         self._Cost = self._results['Cost']
    #         self._Design = self._results['Design']
    #         self._GHGs = self._results['GHG']
    #     else:
    #         raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def _checkbounds(self, key, value, units, bounds):
        """Issue a warning if value is out of bounds.
        
        **Parameters**
        
            **key:** [str] Name of value.
            
            **value:** [float]
            
            **units:** [str] Units of value
            
            **bounds:** iterable[float, float] Upper and lower bounds.
            
        """
        # Warn when value is out of bounds
        lb, ub = bounds
        if not lb<=value and ub>=value:
            units = ' ' + units if units else ''
            try:
                msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb:.4g} to {ub:.4g}{units})."
            except:  # Handle format errors
                msg = f"{key} ({value:.4g}{units}) is out of bounds ({lb} to {ub}{units})."
            warn(_warning(self, msg, DesignWarning), stacklevel=3)
    def _lb_warning(self, key, value, units, lb):
        """Warn that value is below lower bound.
        
         **Parameters**
        
            **key:** [str] Name of value.
            
            **value:** [float]
            
            **units:** [str] Units of value.
            
            **lb:** [float] lower bound.
        
        """
        _lb_warning(key, value, units, lb, 4, self)

    @property
    def ID(self):
        """Unique Identification (str). If set as '', it will choose a default ID."""
        return self._ID

    @ID.setter
    def ID(self, ID):
        if ID == '':
            # Select a default ID if requested
            letter, number = self._default_ID
            self._default_ID[1] += 1
            ID = letter + str(number)
            self._ID = ID
            find.unit[ID] = self
        elif ID and ID != self._ID:
            ID = ID.replace(' ', '_')
            ID_words = ID.split('_')
            if not all(word.isalnum() for word in ID_words):
                raise ValueError('ID cannot have any special characters')
            self._ID = ID
            find.unit[ID] = self

    # Input and output streams
    @property
    def ins(self):
        # list of input streams
        return self._ins    
    @property
    def outs(self):
        # list of output streams
        return self._outs

    @property
    def _downstream_units(self):
        """Return set of all units downstreasm."""
        downstream_units = set()
        outer_periphery = set()
        _add_downstream = _add_downstream_neighbors
        _add_downstream(self, outer_periphery)
        inner_periphery = None
        old_length = -1
        new_length = 0
        while new_length != old_length:
            old_length = new_length
            inner_periphery = outer_periphery
            downstream_units.update(inner_periphery)
            outer_periphery = set()
            for unit in inner_periphery:
                _add_downstream(unit, outer_periphery)
            new_length = len(downstream_units)
        return downstream_units

    def _neighborhood(self, radius=1):
        """Return all neighboring units within given radius.
        
        **Parameters**
        
            **radius:**[int] Maxium number streams between neighbors.
        
        """
        radius -= 1
        neighborhood = set()
        if radius < 0: return neighborhood
        _add_upstream_neighbors(self, neighborhood)
        _add_downstream_neighbors(self, neighborhood)
        direct_neighborhood = neighborhood
        for i in range(radius):
            neighbors = set()
            for neighbor in direct_neighborhood:
                _add_upstream_neighbors(neighbor, neighbors)
                _add_downstream_neighbors(neighbor, neighbors)
            if neighbors == direct_neighborhood: break
            direct_neighborhood = neighbors
            neighborhood.update(direct_neighborhood)
        
        return neighborhood

    def diagram(self, radius=0, file=None, format='png'):
        """Display a `Graphviz <https://pypi.org/project/graphviz/>`__ diagram of the unit and all neighboring units within given radius.
        
        **Parameters**
        
            **radius:** [int] Maxium number streams between neighbors.
        
            **file:** Must be one of the following:
                * [str] File name to save diagram.
                * [None] Display diagram in console.
        
            **format:** Format of file.
        
        """
        if radius > 0:
            neighborhood = self._neighborhood(radius)
            neighborhood.add(self)
            sys = bst.System('', neighborhood)
            return sys.diagram('thorough', file, format)
        
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
            if not stream: continue
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
            if not stream: continue
            f.node(stream.ID) 
            edge_out = self._graphics.edge_out  
            if oi >= len(edge_out): oi = 0
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])
            f.edge(name, stream.ID)
            oi += 1
        save_digraph(f, file, format)
    
    ### Net input and output flows ###
    
    # Molar flow rates
    @property
    def _mol_in(self):
        """Molar flows going in (kmol/hr)."""
        return sum(s.mol for s in self._ins)

    @property
    def _mol_out(self):
        """Molar flows going out (kmol/hr)."""
        return sum(s.mol for s in self._outs)

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
        return sum(s.molnet for s in self._ins)

    @property
    def _molnet_out(self):
        """Net molar flow going out (kmol/hr)."""
        return sum(s.molnet for s in self._outs)

    # Mass flow rates
    @property
    def _mass_in(self):
        """Mass flows going in (kg/hr)."""
        return sum(s.mass for s in self._ins)

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
        return sum(s.massnet for s in self._ins)

    @property
    def _massnet_out(self):
        """Net mass flow going out (kg/hr)."""
        return sum(s.massnet for s in self.outs)

    # Volumetric flow rates
    @property
    def _vol_in(self):
        """Volumetric flows going in (m3/hr)."""
        return sum(s.vol for s in self._ins)

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
        return sum([s.vol for s in self.outs])

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
        return sum([s.H for s in self._ins])

    @property
    def _H_out(self):
        """Enthalpy flow going out (kJ/hr)."""
        return sum([s.H for s in self._outs])

    @property
    def _Hf_in(self):
        """Enthalpy of formation flow going in (kJ/hr)."""
        return sum([s.Hf for s in self._ins])

    @property
    def _Hf_out(self):
        """Enthalpy of formation flow going out (kJ/hr)."""
        return sum([s.Hf for s in self._outs])

    @property
    def _Hnet(self):
        """Net enthalpy flow (including enthalpies of formation)."""
        return self._H_out - self._H_in + self._Hf_out - self._Hf_in
    
    # Representation
    def _info(self, T, P, flow, fraction):
        """Information on unit."""
        if self.ID:
            info = f'{type(self).__name__}: {self.ID}\n'
        else:
            info = f'{type(self).__name__}\n'
        info+= f'ins...\n'
        i = 0
        for stream in self._ins:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            stream_info = stream._info(T, P, flow, fraction)
            unit = stream._source
            index = stream_info.index('\n')
            source_info = f'  from  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + source_info + stream_info[index+1:] + '\n'
            i += 1
        info += f'outs...\n'
        i = 0
        for stream in self._outs:
            if not stream:
                info += f'[{i}] {stream}\n'
                i += 1
                continue
            stream_info = stream._info(T, P, flow, fraction)
            unit = stream._sink
            index = stream_info.index('\n')
            sink_info = f'  to  {type(unit).__name__}-{unit}\n' if unit else '\n'
            info += f'[{i}] {stream.ID}' + sink_info + stream_info[index+1:] + '\n'
            i += 1
        info = info.replace('\n ', '\n    ')
        return info[:-1]

    def show(self, T=None, P=None, flow=None, fraction=None):
        """Prints information on unit."""
        print(self._info(T, P, flow, fraction))
    
    def _ipython_display_(self):
        try: self.diagram()
        except: pass
        self.show()
    
    def _disconnect(self):
        for i in self._ins:
            if i: 
                if i._source: i._sink = None
                else: object.__delattr__(find.stream, i._ID)
        for i in self._outs:
            if i:
                if i._sink: i._source = None
                else: object.__delattr__(find.stream, i._ID)
        self._outs.clear()
        self._ins.clear()
    
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
