# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 11:47:35 2019

@author: yoelr
"""
from numpy import inf

__all__ = ('nearest_neighbor_brute_force',)

def nearest_neighbor_brute_force(x, data):
    diff0 = inf
    sum_ = sum
    abs_ = abs
    for i, y in enumerate(data):
        diffi = sum_(abs_(y - x))
        if diffi < diff0: index = i
    return index