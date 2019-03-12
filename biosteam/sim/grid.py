# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:17:38 2019

@author: yoelr
"""
from biosteam import Unit, np, pd

__all__ = ('Grid',)

DF = pd.DataFrame

# %% Functions

def _blockunit(blockfunc_args):
    element = blockfunc_args[0]._block._element
    if isinstance(element, Unit): return element
    else: return element._sink[0]

def _fill_space(args, argsvalues, argspace, indices):
    for i, test in enumerate(argsvalues):
        for t in test:
            args[i] = t
            last = tuple(args)
            if last not in argspace:
                argspace.append(last)
                indices.append(i) 
                _fill_space(args, argsvalues[:i],
                            argspace, indices)
        test.reverse()


# %% Grid of simulation blocks

class Grid:
    """Create a Grid object that allows for optimized successive simulation of block functions over an argument space.
    
    **Parameters**
    
        **blockfunc_argspace:** [(function, args)] Iterable of block fuctions and respective arguments.
    
    **Examples**

         :doc:`Grid Example`
    
    """
    __slots__ = ('_system',
                 '_funcs',
                 '_initial_args',
                 '_argspace',
                 '_args_array',
                 '_cases',
                 '_name')
    
    options = {'Print exceptions': True}
    
    def __init__(self, system, *blockfunc_argspace, name=None):
        if not blockfunc_argspace: raise ValueError('Must pass at least one blockfunc_argspace')
        
        # Order blocktests from last in network to first
        unitorder = system._flattened_network
        func_args = sorted(blockfunc_argspace,
                            key=lambda x: unitorder.index(_blockunit(x))
                                          if x[0]._block._system else -1,
                            reverse=True)
        
        initial_args = []
        funcs = []
        element_names = []
        argsvalues = []
        for func, args in func_args:
            initial_args.append(args[0])
            funcs.append(func)
            element = func._block._element
            element_names.append(f'{type(element).__name__}-{element}')
            if not isinstance(args, list): args = list(args)
            argsvalues.append(args)
        args = initial_args
        initial_args = tuple(initial_args)
        argspace = []; indices = []
        _fill_space(args, argsvalues, argspace, indices)
        param_names = [i._param.capitalize().replace('_', ' ')
                       for i in funcs]
        cases = [(funcs[i], args[i]) for i, args in
                 zip(indices, argspace)]
        length = len(cases)
        self._args_array = np.asarray(argspace)
        argspace.insert(0, tuple(element_names))
        self._funcs = funcs
        self._system = system
        self._initial_args = initial_args
        if not name:
            self._name = name = funcs[0].__name__
        else:
            self._name = name
        self._cases = cases
        self._argspace = ds = DF(argspace,
                                       columns=param_names,
                                       index=('Unit', *range(length)))
        ds.columns.name = 'Parameter'
        ds[name] = ('n/a', *(None,)*(length))
        
        
    def simulate(self):
        """Simulate all block functions over the argument space."""
        self._argspace[self._name] = ('n/a', *self)
    
    def __iter__(self):
        # Setup units before simulation
        funcs = self._funcs
        for func, case in zip(funcs, self._initial_args):
            func._setter(case)
            e = func._block._element
            if isinstance(e, Unit): e._setup()
        
        # Simulate system to initialize whole network
        self._system.simulate()
        print_exceptions = self.options['Print exceptions']
        for i, (block, case) in enumerate(self._cases):
            try: yield block(case)
            except Exception as e:
                if print_exceptions:
                    d = self._args_array[i, :]
                    e = f'{e}\n'
                    e += str(pd.Series(dict(zip(self._argspace.columns, d))))
                    e = e.rsplit('\n', 2)[0]
                    e += '\n'
                    print(e)
                yield None
    
    @property
    def argspace(self):
        """[DataFrame] Table of the argument space with results in the final column."""
        return self._argspace
        
    def _repr(self):    
        blocks = ''
        for i in self._funcs:
            e = i._block._element
            blk = f'{type(e).__name__}-{e}|'
            if blk not in blocks: blocks += blk
        blocks = blocks.rstrip('|')
        name = self._name
        return f'{type(self).__name__}: {name}'
    
    def __repr__(self):
        return f'[{self._repr()}]'
       
    def _info(self):
        blocks = {}
        for i in self._funcs:
            e = i._block._element
            p = i._param
            blk = f'{type(e).__name__}-{e}'
            if blk not in blocks: blocks[blk] = [p]
            else: blocks[blk].append(p)
        
        lenghts_block = []
        lines = []
        for blk, params in blocks.items():
            newline = f' [{blk}!'
            blklen = len(blk)
            blklen_spaces = (blklen)*' '
            many = len(params) > 1
            if many: newline += '('
            len_params = 0
            for p in params:
                p = f'{p}, '
                len_p = len(p)
                len_params += len_p
                if len_p != len_params and len_params > 30:
                    len_params = 0
                    p = f']\n% [{blklen_spaces}$' + p
                newline += p
            newline = newline.rstrip(', ') + (')' if many else '') + ']\n'
            newlines = newline.split('%')
            lines.extend(newlines)
            lenghts_block.extend((blklen,)*len(newlines))
                
        maxlen_block = max(lenghts_block)
        lengths = []
        lines_2 = []
        for newline, len_ in zip(lines, lenghts_block):
            newline = newline.replace('!', ' '*(maxlen_block-len_) + '| param=')
            newline = newline.replace('$', ' '*(maxlen_block-len_) + '|        ')
            lengths.append(newline.find(']'))
            lines_2.append(newline)
        
        maxlen = max(lengths)
        out = f'{self._repr()}\n'
        for newline, len_ in zip(lines_2, lengths):
            out += newline.replace(']', ' '*(maxlen-len_) + ']')
        return out.rstrip('\n ')
    
    def show(self):
        """Return information on grid."""
        print(self._info())
        
        
        
        
        
        
        
        