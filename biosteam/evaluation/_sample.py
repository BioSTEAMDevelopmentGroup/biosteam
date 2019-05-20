# -*- coding: utf-8 -*-
"""
Created on Sun May 12 17:15:43 2019

@author: yoelr
"""
def cartesian_samples(parameter_values):
    tuple_ = tuple
    initial_sample = []
    sorted_ = sorted
    paramvals = []
    for args in parameter_values:
        args = sorted_(args)
        paramvals.append(args)
        initial_sample.append(args[0])
    num_args = tuple_(enumerate(paramvals))
    sample = initial_sample
    initial_sample = tuple_(initial_sample)
    samples = [initial_sample]
    _fillspace(sample, num_args, samples, tuple_)
    return samples

def _fillspace(args, num_args, argspace, tuple_):
    for i, fargs in num_args:
        for a in fargs[1:]:
            args[i] = a
            argspace.append(tuple_(args))
            _fillspace(args, num_args[:i], argspace, tuple_)
        fargs.reverse()

# def _split(samples, i):
#     subgrids = {}
#     for s in samples:
#         k = s[i]
#         if k in subgrids: subgrids[k].append(s)
#         else: subgrids[k] = [s]
#     return subgrids.values()

# def _sort(grid, i, end):
#     if i == end: return grid
#     subgrids = _split(grid, i)
#     i_next = i+1
#     key = lambda x: x[i_next]
#     grid.clear()
#     reverse = False
#     for samples in subgrids:
#         samples.sort(key=key, reverse=reverse)
#         reverse = not reverse
#         grid.extend(_sort(samples, i_next, end))
#     return grid
    
# def sort(grid):
#     grid.sort(key=lambda x: x[0])
#     end = len(grid[0]) - 1
#     _sort(grid, 0, end)
#     return grid