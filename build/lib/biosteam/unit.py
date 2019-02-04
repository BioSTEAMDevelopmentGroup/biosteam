# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 14:40:28 2018

@author: yoelr
"""

import os
import IPython
from graphviz import Digraph
from biosteam.exceptions import KE
from biosteam.find import find, WeakRefBook
from biosteam.graphics import Graphics, default_graphics
from biosteam.stream import Stream
from biosteam.heat_utility import HeatUtility
from biosteam.utils import get_doc_units, color_scheme, Ins, Outs, MissingStream
from bookkeep import SmartBook, UnitManager
from biosteam.power_utility import PowerUtility
from biosteam import np
CS = color_scheme

dir_path = os.path.dirname(os.path.realpath(__file__)) + '\\'


#%% Decorators and functions

def notify_error(func):
    """Decorate class method to provide a location summary when an error occurs."""
    if hasattr(func, '_original'):
        func = func._original

    def wrapper(self):
        try:
            return func(self)
        except Exception as e:
            # biosteam_Warnings already include location, so it is removed
            location = f'@{type(self).__name__} {self.ID}'
            msg = str(e).strip('\n').replace(location + ': ', '')
            
            # Add location to message
            if not ('@' in msg and ':\n' in msg):
                msg = location + f' {func.__name__}:\n' + msg                 
            
            # Raise exception with same traceback but new message
            import sys
            if type(e) is KeyError:
                raise KE(msg).with_traceback(sys.exc_info()[2])
            else:
                raise type(e)(msg).with_traceback(sys.exc_info()[2])
    
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    wrapper._original = func
    wrapper.__annotations__ = func.__annotations__
    return wrapper


# %% Unit metaclass

Unit_is_done = False
default_line = 'Unit'


class metaUnit(type):
    """Unit metaclass for wrapping up methods with error notifiers, adding key word arguments, and keeping track for Unit lines and inheritance. Also adds the instances attribute to keep track of instances of the class.

    **Class definitions**

        kwargs = {}: [dict] Default keyword arguments

        setup(): Create components and cached data during initialization

        run(): Run rigorous simulation

        operation(): Find operation requirements

        design(): Find design requirements

        cost(): Find capital and annual cost

        .. Note::

           All definitions are optional

    **Class Attributes**
    
        CEPCI = 567.5: Chemical plant cost index

    **Instance Attributes**

        instances = {}: Dictionary of weak references to instances of this class

        line = [Defaults as the class name of the first child class]: [str] Name denoting the type of Unit class

    """
    _CEPCI = 567.5 # Chemical engineering plant cost index (567.5 at 2017)
    def __new__(mcl, clsname, superclasses, new_definitions):
        """Prepare unit methods with wrappers for error notification, and add kwargs as key word arguments to __init__. Also initiallize by adding new Unit class to the find.line dictionary."""

        if not Unit_is_done:
            # Abstract Unit class
            cls = type.__new__(mcl, clsname, superclasses, new_definitions)
        else:
            # Add error notification wrapper to every new unit method
            for method_name in ('setup', 'run', 'operation', 'design', 'cost'):
                if new_definitions.get(method_name):
                    new_definitions[method_name] = notify_error(new_definitions[method_name])
    
            # Make new Unit class
            cls = type.__new__(mcl, clsname, superclasses, new_definitions)
    
            if cls.__doc__ is Unit.__doc__:
                # Do not inherit docstring from Unit
                cls.__doc__ = None
            if cls.__doc__:
                __doc__ = "*An extension of the* :doc:`Unit class <Unit>`. " + cls.__doc__
                cls.__doc__ = __doc__.replace('**Parameters**', 
            "**Parameters**\n\n" +
            
        "        **ID:** [str] Unique identification. If set as '', a default ID will be chosen.\n\n" +

        "        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.\n\n" +
        
        "        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.")
    
            # Set line
            # default_line constitutes a new Unit class
            if cls.line is default_line:
                cls.line = cls.__name__
                # Set new graphics object for new line
                if not new_definitions.get('_Graphics'):
                    cls._Graphics = Graphics()
            elif new_definitions.get('line'):
                # Set new graphics for specified line
                if not new_definitions.get('_Graphics'):
                    cls._Graphics = Graphics()
    
            # Add class to find line dictionary
            if cls.line is not default_line: # Do not include abstract Unit classes
                if not find.line[cls.line]:
                    find.line[cls.line] = []
                if cls.__name__ not in [cls.__name__ for cls in find.line[cls.line]]:
                    find.line[cls.line].append(cls)
                    
            # Key word arguments to replace
            kwargs = ''
            inputs = ''
            for key, val in cls.kwargs.items():
                if type(val) is str:
                    val = f"'{val}'"
                kwargs += ', ' + key + ' = ' + str(val)
                inputs += ', ' + key + ' = ' + key
    
            # Begin changing __init__ to have kwargs by making a string to execute
            str2exec = f"def __init__(self, ID='', outs=(), ins=None{kwargs}):\n"
            str2exec += f"     initfunc(self, ID, outs, ins{inputs})"
    
            # Execute string and replace __init__
            globs = {'initfunc': Unit.__init__,
                     'Stream': Stream}
            locs = {}
            exec(str2exec, globs, locs)
            cls.__init__ = locs['__init__']

        #: Dictionary of weak references to instances of this class
        cls._instances = WeakRefBook(check_ID=False)
        for method_name in ('operation', 'design', 'cost'):
            method = getattr(cls, method_name)
            annotations = method.__annotations__
            if not annotations:
                annotations['return'] = method_name.capitalize() + ' [dict]'

        # Get units of measure of each key method
        units = {}
        get_doc_units(units, cls.operation.__doc__)
        get_doc_units(units, cls.design.__doc__)
        get_doc_units(units, cls.cost.__doc__)
        cls._results_UnitsOfMeasure = UnitManager([], **units)
        
        return cls
    
    @property
    def CEPCI(cls):
        """Chemical engineering plant cost index (CEPCI)"""
        return cls._CEPCI
    @CEPCI.setter
    def CEPCI(cls, CEPCI):
        cls._CEPCI = CEPCI
    
    @property
    def instances(cls):
        """Dictionary of weak references to all instances of the Unit class."""
        return cls._instances
    
    def __repr__(cls):
        if cls.line == 'Unit' and cls.__name__ == 'Unit':
            return f'biosteam.{cls.__name__}'
        elif cls.line == cls.__name__ or cls.line == 'Unit':
            return f'Unit.{cls.__name__}'
        else:
            return f'{cls.line}.{cls.__name__}'

# %% Unit Operation


class Unit(metaclass=metaUnit):
    """Abstract parent class for Unit objects. Child objects must contain setup, run, operation, design and cost methods to setup internal objects, estimate stream outputs of a Unit and find operating, design and cost information. These methods should store information in the 'results' dictionary (an instance attribute).  

    **Parameters**

        **ID:** [str] Unique identification

        **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.
        
        **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.

        ****kwargs:** Keyword arguments that are accessed by setup, run, operation, design, and cost methods

    **Abstract Attributes**

        **kwargs** = {}: Default keyword arguments.

        **setup()**
            Create components and cached data. 

        **run()**
            Run rigorous simulation.

        **operation()** -> Operation [dict]
            Update results and return dictionary of operation requirements.

        **design()** -> Design [dict]
            Update results and return dictionary of design requirements.

        **cost()** -> Cost [dict]
            Update results and return dictionary of costs.

    **Class Attributes** 

        **CEPCI** = 567.5: [float] Chemical Engineering Plant Cost Index

        **operating_days**  = 330: [int] Operation days per year
        
        **lang_factor** = 5.03: [float] Lang factor including working capital from Peters, Timmerhaus, and West (2003)
        
        **bounds** = {} [dict] Values should be tuples with lower and upper bounds.
        
        **instances** = {}: [WeakRefBook] Contains all instances of the Unit class.
        
        **_N_ins** = 1: [int or None] Expected number of input streams

        **_N_outs** = 2: [int or None] Expected number of output streams

        **_N_heat_util** = 0: [int] Number of heat utilities  

        **_power_util** = False: [bool] If True, a PowerUtility object is created for every instance.
        

        .. Note:

           The Unit class is an instance of metaUnit, which takes class definitions (defined in abstract attributes) to set key word arguments, decorate methods, and register new classes.

    **ins**
        
        list of input streams
        
    **outs**
    
        list of output streams
    
    """ 
    ### Abstract Attributes ###
    
    # [int or None] Expected number of input streams
    _N_ins = 1  
    
    # [int or None] Expected number of output streams
    _N_outs = 2  
    
    # [dict] Values should be tuples with lower and upper bounds for results dictionary.
    bounds = {}
    
    # [float] Lang factor including working capital from Peters, Timmerhaus, and West (2003)
    lang_factor = 5.03 
    
    # [int] Operation days per year
    operating_days = 330
    
    # [int] number of heat utilities
    _N_heat_util = 0
    
    # [PowerUtility] A PowerUtility object, if required
    _power_utility = None
    
    # [bool] If True, a PowerUtility object is created for every instance.
    _power_util = False 
    
    # [dict] default key word arguments that are accessed by setup, run, _simple_run, operation, design, and cost methods.
    kwargs = {}
    
    # [biosteam Graphics] a Graphics object for diagram representation.
    _Graphics = default_graphics
    
    # [str] The type of unit, regardless of estimations
    line = default_line

    ### Other defaults ###
    
    # [float] Contingency and fees factor :math:`C_{TM} = F_{TM} * C_{BM}`
    _F_TM = 1.18

    # [list] Default ID starting letter and number
    _default_ID = ['U', 0]
    
    # Default ID
    _ID = ''  
    
    #: [bool] True if unit is in a recycle loop
    _in_loop = False

    ### Initialize ###
    
    def __init__(self, ID='', outs=(), ins=None, **kwargs):
        """Initialize Unit object. See help(type(self)) for accurate signature.

        **Parameters**

            **ID:** [str] Unique identification

            **outs:** tuple[str or Stream] Output streams or IDs to initialize output streams. If None, leave streams missing. If empty, default IDs will be given.
            
            **ins:** tuple[str or Stream] Input streams or IDs to initialize input streams. If None, leave streams missing. If empty, default IDs will be given.
    
            ****kwargs:** Keyword arguments that are accessed by setup, run, operation, design, and cost methods
        
        """
        self.ID = ID
        self.kwargs = kwargs  #: [dict] Key word arguments
        self._init_ins(ins)
        self._init_outs(outs)
        self._init_results()
        self._init_heat_utils()
        self._init_power_util()
        self.setup()

    def _init_ins(self, ins):
        # Initialize input streams
        if ins is None:
            self._ins = Ins(self, (MissingStream() for i in range(self._N_ins)))
        elif isinstance(ins, Stream):
            self._ins = Ins(self, (ins))
        elif isinstance(ins, str):
            self._ins = Ins(self, (Stream(ins)))
        elif not ins:
            self._ins = Ins(self, (Stream('') for i in range(self._N_ins)))
        else:
            self._ins = Ins(self, (Stream(i) if isinstance(i, str) else i for i in ins))
    
    def _init_outs(self, outs):
        # Initialize output streams
        if isinstance(outs, str):
            self._outs = Outs(self, (Stream(outs)))
        elif isinstance(outs, Stream):
            self._outs = Outs(self, (outs))
        elif outs is None:
            self._outs = Outs(self, (MissingStream() for i in range(self._N_outs)))
        elif not outs:
            self._outs = Outs(self, (Stream('') for i in range(self._N_outs)))
        else:
            self._outs = Outs(self, (Stream(i) if isinstance(i, str) else i for i in outs))
    
    def _init_results(self):
        # Initialize results attribute
        units = self._results_UnitsOfMeasure
        bounds = self.bounds
        empty = {}
        self._results = SmartBook(units, bounds, source=self, inclusive=empty,
            Operation=SmartBook(units, bounds, source=self, inclusive=empty),
            Design=SmartBook(units, bounds, source=self, inclusive=empty),
            Cost=SmartBook(units, bounds, source=self, inclusive=empty))
    
    def _init_heat_utils(self):
        # Initialize heat utilities
        _N_heat_util = self._N_heat_util
        if _N_heat_util:
            utilities = [HeatUtility() for i in range(_N_heat_util)]
        else:
            utilities = None
        self._heat_utilities = utilities
        
    def _init_power_util(self):
        # Initialize power utility
        if self._power_util:
            self._power_utility = PowerUtility()
    
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
    def __neg__(self):
        return self.ins
    
    # Abstract methods
    def setup(self): pass
    def _run(self) pass
    def operation(self)->'Operation [dict]': return self.results['Operation']
    def design(self)->'Design [dict]': return self.results['Design']
    def cost(self)->'Cost [dict]': return self.results['Cost']

    def simulate(self):
        """Run simulation, operation, design and cost methods."""
        self.run()
        self.operation()
        self.design()
        self.cost()

    @property
    def results(self):
        """[dict] Key results from running Unit methods."""
        return self._results

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
        # Remove old reference to this object
        if self._ID != '' and ID != '':
            del find.unit[self._ID]
            del type(self).instances[self._ID]

        # Get current default ID
        Unit = type(self)
        letter, number = Unit._default_ID

        # Make sure given ID is not a default ID
        if ID.startswith(letter):
            if ID[1:].isdigit():
                raise ValueError(f"IDs starting with '{letter}' and followed by only digits are reserved for defaults. To use the default ID, set ID to 'Default'.")

        # Select a default ID if requested
        if ID == '':
            Unit._default_ID[1] += 1
            ID = letter + str(number)

        # Add instance to class instance dictionary
        type(self).instances[ID] = self

        # Add ID to find dictionary
        find.unit[ID] = self
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

    # Utilities
    @property
    def heat_utilities(self):
        """list of HeatUtility objects"""
        return self._heat_utilities
    @heat_utilities.setter
    def heat_utilities(self, heat_utilities):
        self._heat_utilities = heat_utilities

    @property
    def power_utility(self):
        """A PowerUtility object, if required."""
        return self._power_utility
    @power_utility.setter
    def power_utility(self, power_utility):
        if not isinstance(power_utility, PowerUtility):
            raise TypeError(f'Only PowerUtility objects are valid.')
        else:
            self._power_utility = power_utility

    @property
    def diagram(self):
        """`Graphviz diagram <https://pypi.org/project/graphviz/>`__ of unit"""

        # Make a Digraph handle
        f = Digraph(name='unit', filename='unit', format='svg')
        f.attr('graph', ratio='0.5', splines='normal', outputorder='edgesfirst',
               nodesep='1.1', ranksep='0.8', maxiter='1000')  # Specifications
        f.attr(rankdir='LR')  # Left to right

        # If many streams, keep streams close
        if (len(self.ins) >= 3) or (len(self.outs) >= 3):
            f.attr('graph', nodesep='0.4')

        # Names for unit nodes
        type_ = type(self).__name__.replace('_', '\n').replace(' ', '\n')
        name = self.ID + '\n' + type_

        # Initialize node arguments based on unit and make node
        self._Graphics.node_function(self)
        f.attr('node', **self._Graphics.node)
        f.node(name)

        # Set stream node attributes
        f.attr('node', shape='rarrow', fillcolor='#79dae8',
               style='filled', orientation='0', width='0.6',
               height='0.6', color='black')

        # Make nodes and edges for input streams
        di = 0  # Destination position of stream
        for stream in self.ins:
            # If stream is not important (no ID), ignore the stream
            if stream.ID == '' or stream.ID == 'Missing Stream':
                continue
            f.node(stream.ID)  # Make node
            edge_in = self._Graphics.edge_in  # Get edge attributes
            f.attr('edge', arrowtail='none', arrowhead='none',
                   tailport='e', **edge_in[di])  # Set edge attibutes
            f.edge(stream.ID, name)  # Make edge
            di += 1

        # Make nodes and edges for output streams
        oi = 0  # Origin position of stream
        for stream in self.outs:
            f.node(stream.ID)  # Make node
            edge_out = self._Graphics.edge_out  # Get edge attributes
            f.attr('edge', arrowtail='none', arrowhead='none',
                   headport='w', **edge_out[oi])  # Set edge attibutes
            f.edge(name, stream.ID)  # Make edge
            oi += 1

        # Display digraph on console
        x = IPython.display.SVG(f.pipe(format='svg'))
        IPython.display.display(x)

    # Costs
    @property
    def operating_cost(self) -> 'USD/hr':                
        """Sum of all utility costs."""
        cost = 0
        heat_utilities = self.heat_utilities
        power_utility = self.power_utility
        
        if heat_utilities:
            cost = sum(util.results['Cost'] for util in self.heat_utilities)
        if power_utility:
            cost += self.power_utility.results['Cost']
        return cost
        
    @property
    def capital_cost(self) -> 'USD':
        """The sum of all costs calculated by the cost method."""
        return sum(self.results['Cost'].values())
    
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
            # Get non-zero index and species
            index = []
            species = []
            for i in stream._index:
                if stream.mol[i] != 0:
                    index.append(i)
                    species.append(stream._index_ID[i])
            # Find Hf of stream
            Hf = sum(stream.mol[index] *
                     stream._species.get_props(species, 'Hfm'))
            _Hf_in += Hf
        return _Hf_in

    @property
    def _Hf_out(self):
        """Enthalpy of formation flow going out (kJ/hr)."""
        _Hf_out = 0
        for stream in self.outs:
            # Get non-zero index and species
            index = []
            species = []
            for i in stream._index:
                if stream.mol[i] != 0:
                    index.append(i)
                    species.append(stream._index_ID[i])
            # Find Hf of stream
            Hf = sum(stream.mol[index] *
                     stream._species.get_props(species, 'Hfm'))
            _Hf_out += Hf
        return _Hf_out

    @property
    def Hnet(self):
        """Net enthalpy flow (including enthalpies of formation)."""
        return self._H_out - self._H_in + self._Hf_out - self._Hf_in
    
    # Representation
    def _info(self, **show_units):
        """Information on unit."""
        info = (f'{type(self).__name__}: {self.ID}\n'
                + f'{CS.dim("ins...")}\n')
        i = 0
        for stream in self.ins:
            if stream.ID == 'Missing Stream':
                info += f'[{i}] {stream.ID}\n'
                i += 1
                continue
            unit = stream._source[0]
            stream_info = stream._info(**show_units)
            index = stream_info.index('\n')
            if unit is None:
                source_info = '\n'
            else:
                source_info = f'  from  {unit}\n'
            info += f'[{i}] {stream.ID}' + source_info + stream_info[index+1:] + '\n'
            i += 1
        info += f'{CS.dim("outs...")}\n'
        i = 0
        for stream in self.outs:
            if stream.ID == 'Missing Stream':
                info += f'[{i}] {stream.ID}\n'
                i += 1
                continue
            unit = stream._sink[0]
            stream_info = stream._info(**show_units)
            index = stream_info.index('\n')
            if unit is None:
                sink_info = '\n'
            else:
                sink_info = f'  to  {unit}\n'
            info += f'[{i}] {stream.ID}' + sink_info + stream_info[index+1:] + '\n'
            i += 1
        info = info.replace('\n ', '\n    ')
        return info[:-1]

    def show(self, **show_units):
        """Prints information on unit."""
        print(self._info(**show_units))
    
    def __str__(self):
        return self.ID

    def __repr__(self):
        return '<' + type(self).__name__ + ': ' + self.ID + '>'

Unit_is_done = True
