# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
import numpy as np
import pandas as pd
from .. import Unit
from ._block import Block

__all__ = ('Grid',)

DF = pd.DataFrame

# %% Functions

def _blockunit(blockfunc_args):
    element = blockfunc_args[0]._block._element
    if isinstance(element, Unit): return element
    else: return element._sink

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

class Grid:
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
                 '_ID',
                 '_layout')
    
    def __init__(self, system, metric, *, args=(), ID=None, layout='simple'):
        if args:
            metric_no_args = lambda: metric(*args)
            metric_no_args.__name__ = metric.__name__
            self._metric = metric_no_args
        else: self._metric = metric
        self.layout = layout
        self._system = system
        self._ID = ID or self._metric.__name__
        
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
    
    def evalparam(self, element, setter, values, isolated=False):
        """Return metric at given parameter values.
        
        **Parameters**
        
            **element:** [Unit or Stream] Element in the system being altered.
            
            **setter:** [function] Should set parameter in the element.
            
            **values:** [iterable] Values for parameter.
            
            **isolated** [bool] If True, account for downstream operations. If False, only account for element.
        
        """
        system = None if isolated else self._system
        blockfunc = Block(element, system)(self._metric, setter)
        return [blockfunc(i) for i in values]
    
    def addparam(self, element, setter, values, isolated=False):
        """Add parameter to vary in metric.
        
        **Parameters**
        
            **element:** [Unit or Stream] Element in the system being altered.
            
            **setter:** [function] Should set parameter in the element.
            
            **values:** [iterable] Values for parameter.
            
            **isolated** [bool] If True, account for downstream operations. If False, only account for element.
        
        """
        system = None if isolated else self._system
        self._stack.append((Block(element, system)(self._metric, setter), values))
    
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
        if not stack: raise RuntimeError(f'No metric parameters set for {repr(self)} object.')
        tuple_ = tuple
        # Order blocktests from last in network to first
        length = len(system._unitnetwork)
        index = system._unitnetwork.index
        stack.sort(key=lambda x: index(_blockunit(x))
                                 if x[0]._block._system else length,
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
                element = func._block._element
                element_names.append(element.line + '-' + element.ID.replace('_', ' '))
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
                element = func._block._element
                element_names.append(element.line + '-' + element.ID.replace('_', ' '))
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
        
    def simulate(self):
        """Simulate Grid object over the argument space and save metric values to `table`."""
        # Setup units before simulation
        argspace, funcs, initial_args, thread = self._loadtable()
        for func, arg in zip(funcs, initial_args):
            func._setter(arg)
        
        # Simulate system to initialize whole network
        try: self._system.simulate()
        except: pass
        values = []
        add = values.append
        if thread:
            for func, arg in thread:
                try: add(func(arg))
                except: add(None)
        else:
            metric = self._metric
            add(metric())
            last_args = initial_args
            simulate_funcs = [i._block._simulate for i in funcs]
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
        
    def _repr(self):    
        return f'{type(self).__name__}: {self._ID}'
    
    def __repr__(self):
        return f'<{self._repr()}>'
       
    def _info(self):
        if not self._stack:
            return (f'{self._repr()}\n'
                    +' element  parameters\n'
                    +' None     None')
        blocks = {}
        for i, _ in self._stack:
            e = i._block._element
            p = i._param
            blk = f'{type(e).__name__}-' + f'{e}'.replace('_', ' ')
            if blk not in blocks: blocks[blk] = [p]
            else: blocks[blk].append(p)
        
        lenghts_block = []
        lines = []
        for blk, params in blocks.items():
            blklen = len(blk)
            blklen_spaces = (blklen)*' '
            newlines = []
            newlines.append(f" {blk}${params[0].replace('_', ' ')}\n")
            for p in params[1:]:
                newlines.append(f" {blklen_spaces}${p.replace('_', ' ')}\n")
            lines.extend(newlines)
            lenghts_block.extend((blklen,)*len(newlines))
        
        maxlen_block = max(lenghts_block)
        out = f'{self._repr()}\n'
        maxlen_block = max(maxlen_block, 7)
        out += ' element' + (maxlen_block - 7)*' ' + '  parameter\n'
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '  ')
            out += newline
        
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on metric parameters."""
        print(self._info())
        
        
        