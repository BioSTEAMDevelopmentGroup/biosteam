# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
import numpy as np
import pandas as pd
from ._model import elementname, blockunit, blockfunction, Model

__all__ = ('Grid',)

DF = pd.DataFrame

# %% Functions

def _fillspace(args, num_stack, argspace, thread, tuple_):
    for i, f, fargs in num_stack:
        for a in fargs[1:]:
            args[i] = a
            argspace.append(tuple_(args))
            thread.append((f, a))
            _fillspace(args, num_stack[:i], argspace, thread, tuple_)
        fargs.reverse()

def _orderedspace(args, num_args, argspace, tuple_, space):
    for i, fargs in num_args:
        for a in fargs[1:]:
            args[i] = a
            new = tuple_(args)
            if new in space:
                argspace.append(new)
            _orderedspace(args, num_args[:i], argspace, tuple_, space)
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
    __slots__ = ('_ID',     # [str] Should be the metric name.
                 '_table',  # [DataFrame] All arguments and results.
                 '_system', # [System] Reflects the model state.
                 '_stack',  # list[(function, args)] Model block fuctions and respective arguments.
                 '_params', # [list] Cached parameter values of last model run.
                 '_metric', # [function] Should return metric being evaluated.
                 '_blockf', # list[function] All block functions in the stack.
                 '_cached', # tuple[argspace, blockf, initial_args, thread] Cached values from loading table.
                 '_reload', # [bool] True if table should be reloaded.
                 '_layout') # [str] {'simple', 'cartesian'}
    
    def __init__(self, ID, system, metric, *, args=(), layout='simple'):
        super().__init__(system)
        self._ID = ID
        self._metric = (lambda: metric(*args)) if args else metric
        self.layout = layout
        #: list[(function, args)] Block fuctions and respective arguments.
        self._stack = []
        #: [DataFrame] Table of the argument space with results in the final column.
        self._table = None
        self._cached = 4*(None,)
    
    def __len__(self):
        return len(self._stack)
    
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
        if self._reload: self._loadtable()
        return self._table
    
    def addparam(self, setter=None, element=None, values=None, kind='isolated', param=None):
        """Add parameter to vary in metric.
        
        **Parameters**
        
            **setter:** [function] Should set parameter in the element.
            
            **element:** [Unit or Stream] Element in the system being altered.
            
            **values:** [iterable] Values for parameter.
            
            **kind:** {'isolated', 'design', 'cost'}
                * 'coupled': parameter is coupled to the system.
                * 'isolated': parameter does not affect the system.
                * 'design': parameter only affects design and cost of the element.
                * 'cost': parameter only affects cost of the element.
                
            **param:** [str] Name of parameter. If None, default to argument name of setter.
            
        .. Note::
            
            If kind is 'coupled', account for downstream operations. Otherwise, only account for given element. If kind is 'design' or 'cost', element must be a Unit object.
        
        """
        if not setter: return lambda setter: self.addparam(setter, element, values, kind, param)
        self._stack.append((blockfunction(self._system, element, setter, kind, param), values))
        self._reload = True

    def _loadmodel(self):
        if self._reload:    
            # Order blocktests from last in network to first
            length = len(self._system._unitnetwork)
            index = self._system._unitnetwork.index
            self._stack.sort(key=lambda x: index(blockunit(x[0]))
                                           if x[0]._system else length)
            self._blockf = [i[0] for i in self._stack]
    _loadmodel.__doc__ = Model._loadmodel.__doc__
    
    def _loadtable(self):
        """Load argument space and return parameters for simulation."""
        if not self._reload: return
        self._loadmodel()
        blockf = self._blockf
        tuple_ = tuple
        initial_args = []
        element_names = []
        all_args = []
        list_ = list
        if self._layout == 'cartesian':
            for f, args in self._stack:
                initial_args.append(args[0])
                element_names.append(elementname(f._element))
                all_args.append(list_(args))
            num_stack = tuple_(zip(range(len(blockf)), blockf, all_args))
            args = initial_args
            initial_args = tuple_(initial_args)
            argspace = [initial_args]
            thread = [(blockf[0], initial_args[0])]
            _fillspace(args, num_stack, argspace, thread, tuple_)
        else:
            argspace = []
            for f, args in self._stack:
                argspace.append(args)
                element_names.append(elementname(f._element))
            argspace = list_(zip(*argspace))
            key = lambda x: x[i]
            for i in range(len(blockf)-1,  -1, -1): argspace.sort(key=key)
            args = initial_args = argspace[0]
            thread = None
        argspace = np.asarray(argspace)
        paramIDs = [i._param.capitalize().replace('_', ' ')
                    for i in blockf]
        spacelen = len(argspace)
        spacerange = range(spacelen)
        self._table = DF(argspace,
                         columns=pd.MultiIndex.from_arrays((element_names, paramIDs),
                                                           names=('Element', 'Parameter')),
                         index=spacerange)
        self._table[self._ID] = (None,)*spacelen
        self._cached = (argspace, blockf, initial_args, thread)
        self._reload = False
        
    def evaluate(self, default=None):
        """Evaluate Grid object over the argument space and save metric values to `table`."""
        # Setup before simulation
        self._loadtable()
        argspace, modelf, initial_args, thread = self._cached
        for f, arg in zip(modelf, initial_args): f._setter(arg)
        metric = self._metric
        values = []
        add = values.append
        if thread:
            for func, arg in thread:
                try: func(arg)
                except: add(default)
                finally: add(metric())
        else:
            for args in argspace: add(self(args, default))
        self._table[self._ID] = values
    
    def __call__(self, parameters, default=None):
        super().__call__(parameters, default)
        return self._metric()
    
    def __repr__(self):
        return super().__repr__()[:-1] + f' -> {self._ID}>'
    
    def _repr(self):
        return f'{type(self).__name__}: {self._ID}'
        
   
    
        