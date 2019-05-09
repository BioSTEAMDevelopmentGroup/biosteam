# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
import numpy as np
import pandas as pd
from ._block import Block
from ._model import _get_element_name, _blockunit, Model

__all__ = ('Grid',)

DF = pd.DataFrame

# %% Functions

def _fill_space(args, num_stack, argspace, thread, tuple_):
    for i, f, fargs in num_stack:
        for a in fargs[1:]:
            args[i] = a
            argspace.append(tuple_(args))
            thread.append((f, a))
            _fill_space(args, num_stack[:i], argspace, thread, tuple_)
        fargs.reverse()

def _ordered_space(args, num_args, argspace, tuple_, space):
    for i, fargs in num_args:
        for a in fargs[1:]:
            args[i] = a
            new = tuple_(args)
            if new in space:
                argspace.append(new)
            _ordered_space(args, num_args[:i], argspace, tuple_, space)
        fargs.reverse()

# %% Grid of simulation blocks

class Grid(Model):
    """Create a Grid object that allows for optimized evaluation over an argument space.
    
    **Parameters**
    
        **system:** [System] Metric parameters should act on the system.
    
        **metric:** [function] Should return metric valye.
        
        **args:** [tuple] Arguments to pass to the `metric` function.
        
        **ID:** [str] ID of the Grid object
        
        **layout:** {'simple', 'cartesian'}
            * 'simple': The argument space is a matrix of the parameter values.
            * 'cartesian': The argument space is a cartesian product of the parameter values.
    
    **Examples**

         :doc:`Grid Example`
    
    """
    __slots__ = ('_table',
                 '_system',
                 '_stack',
                 '_metric',
                 '_params',
                 '_orderf',
                 '_ID',
                 '_layout')
    
    def __init__(self, ID, system, metric, *, args=(), layout='simple'):
        super().__init__(ID, system, metric, args=args)
        self.layout = layout
        #: [(function, args)] Iterable of block fuctions and respective arguments.
        self._stack = []
        #: [DataFrame] Table of the argument space with results in the final column.
        self._table = None
    
    @property
    def layout(self):
        """
        Must be either 'simple' or 'cartesian'
            * 'simple': The argument space is a matrix of the parameter values.
            * 'cartesian': The argument space is the cartesian product of parameter values.
        """
        return self._layout
    
    @layout.setter
    def layout(self, layout):
        layout = layout.casefold()
        if layout not in ('simple', 'cartesian'):
            if isinstance(layout, str): layout = f"'{layout}'"
            raise ValueError(f"Grid layout must be either 'simple' or 'cartesian', not {layout}")
        self._layout = layout
    
    @property
    def table(self):
        """[DataFrame] Table of the argument space with results in the final column."""
        if not self._stack: return None
        elif self._table is None: self._loadtable()
        return self._table
    
    def evalparam(self, element, setter, values, kind='isolated'):
        """Return metric at given parameter values.
        
        **Parameters**
        
            **element:** [Unit or Stream] Element in the system being altered.
            
            **setter:** [function] Should set parameter in the element.
            
            **values:** [iterable] Values for parameter.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system.
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if kind is 'coupled':
            blockf = Block(element, self._system)(self._metric, setter)
        elif kind is 'isolated':
            blockf = Block(element, None)(self._metric, setter)
        elif kind is 'design':
            blockf = Block(element, None)(self._metric, setter, element._summary)
        elif kind is 'cost':
            blockf = Block(element, None)(self._metric, setter, element._finalize)
        else:
            raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")
        return [blockf(i) for i in values]
    
    def addparam(self, element, setter, values, kind='isolated'):
        """Add parameter to vary in metric.
        
        **Parameters**
        
            **element:** [Unit or Stream] Element in the system being altered.
            
            **setter:** [function] Should set parameter in the element.
            
            **values:** [iterable] Values for parameter.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system.
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if kind is 'coupled':
            blockf = Block(element, self._system)(setter)
        elif kind is 'isolated':
            blockf = Block(element, None)(setter)
        elif kind is 'design':
            blockf = Block(element, None)(setter, element._summary)
        elif kind is 'cost':
            blockf = Block(element, None)(setter, element._finalize)
        else:
            raise ValueError(f"kind must be either 'coupled', 'isolated', 'design', or 'cost' (not {kind}).")
        self._stack.append((blockf, values))

    @property
    def _blockf(self):
        return [i[0] for i in self._stack]
    @_blockf.setter
    def _blockf(self, value): pass
    
    def _loadtable(self):
        """Load argument space and return parameters for simulation.
        
        **Returns**
        
            **argspace:** [array] All arguments
            **funcs:** [list] All block functions
            **initial_args:** [tuple] First arguments to be simulated
            **thread:** [list] All block functions and arguments to be simulated.
        
        """
        system = self._system
        stack = self._stack
        tuple_ = tuple
        # Order blocktests from last in network to first
        length = len(system._unitnetwork)
        index = system._unitnetwork.index
        stack.sort(key=lambda x: index(_blockunit(x[0]))
                                 if x[0]._system else length,
                   reverse=True)
        
        initial_args = []
        funcs = []
        element_names = []
        all_args = []
        list_ = list
        if self._layout == 'cartesian':
            for func, args in stack:
                initial_args.append(args[0])
                funcs.append(func)
                element = func._element
                element_names.append(_get_element_name(element))
                all_args.append(list_(args))
            num_stack = tuple_(zip(range(len(funcs)), funcs, all_args))
            args = initial_args
            initial_args = tuple_(initial_args)
            argspace = [initial_args]
            thread = [(funcs[0], initial_args[0])]
            _fill_space(args, num_stack, argspace, thread, tuple_)
        else:
            argspace = []
            for func, args in stack:
                argspace.append(args)
                funcs.append(func)
                element = func._element
                e_name = _get_element_name(element)
                element_names.append(e_name)
            argspace = list_(zip(*argspace))
            key = lambda x: x[i]
            index = list_(range(len(stack)))
            index.reverse()
            for i in index: argspace.sort(key=key)
            args = initial_args = argspace[0]
            thread = None
        
        paramIDs = [i._param.capitalize().replace('_', ' ')
                    for i in funcs]
        spacelen = len(argspace)
        spacerange = range(spacelen)
        self._table = DF(argspace,
                         columns=pd.MultiIndex.from_arrays((element_names, paramIDs),
                                                           names=('Element', 'Parameter')),
                         index=spacerange)
        self._table[self._ID] = (None,)*spacelen
        
        return (np.asarray(argspace),
                funcs,
                initial_args,
                thread)
        
    def evaluate(self):
        """Evaluate Grid object over the argument space and save metric values to `table`."""
        # Setup units before simulation
        argspace, funcs, initial_args, thread = self._loadtable()
        for func, arg in zip(funcs, initial_args):
            func._setter(arg)
        
        metric = self._metric
        # Simulate system to initialize whole network
        try: self._system.simulate()
        except: pass
        values = []
        add = values.append
        if thread:
            for func, arg in thread:
                try: 
                    func(arg)
                    add(metric())
                except: add(None)
        else:
            metric = self._metric
            add(metric())
            last_args = initial_args
            simulate_funcs = [i._simulate for i in funcs]
            setter_funcs = [i._setter for i in funcs]
            index = tuple(range(len(last_args)))
            for args in argspace[1:]:
                simfunc_pos = None
                for i, j, k in zip(index, args, last_args):
                    if j != k:
                        setter_funcs[i](j)
                        simfunc_pos = i
                if simfunc_pos: simulate_funcs[simfunc_pos]()
                add(metric())
                last_args = args
        self._table[self._ID] = values
        
    __repr__ = Model.__repr__
    _info = Model._info
    _repr = Model._repr
    show = Model.show
        
        
        