# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module includes arbitrary classes and functions.

"""
import biosteam as bst
from thermosteam import Chemicals
from math import ceil

__all__ = ('factor', 'checkbounds', 'strtuple',
           'format_title', 'format_unit_name',
           'remove_undefined_chemicals',
           'default_chemical_dict', 'subgroup',
           'repr_subgroups', 'repr_items')

# %% Number functions

def factor(base_units, new_units):
    if base_units == new_units: return 1
    else: return bst._Q(1, base_units).to(new_units).magnitude

def checkbounds(x, bounds):
    return bounds[0] < x < bounds[1]

# %% String functions

def strtuple(iterable):
    """Return string of all items in the tuple""" 
    string = ''
    function = type(strtuple)
    for i in iterable:
        if isinstance(i , function):
            string += i.__name__ + ', '
        else:
            string += str(i) + ', '
    string = string.rstrip(', ')
    string = '(' + string + ')'
    return string
        
def format_title(line):
    line = line.replace('_', ' ')
    words = []
    word = ''
    last = ''
    for i in line:
        if i.isupper() and last.isalpha() and not last.isupper():
            words.append(word)
            word = i
        else:
            word += i
        last = i
    words.append(word)
    line = ''
    for word in words:
        N_letters = len(word)
        if N_letters > 1:
            line += word + ' '
        else:
            line += word
    line = line.strip(' ')
    first_word, *rest = line.split(' ')
    words = [first_word[0].capitalize() + first_word[1:]]
    for word in rest:
        if not all([(i.isupper() or not last.isalpha()) for i in word]):
            word = word.lower()
        words.append(word)
    return ' '.join(words)

def format_unit_name(name):
    return ''.join([i[0].capitalize() + i[1:] for i in name.split(' ')])
    
def subgroup(items, size=5):
    return [items[size*i: size*(i+1)] for i in range(int(ceil(len(items) / size)))]

def repr_subgroups(subgroups):
    return [', '.join([str(i) for i in j]) for j in subgroups]

def repr_items(start, items, subgroup_size=5, brackets=None):
    N_spaces = len(start)
    if brackets:
        left, right = brackets
        N_spaces += 1
    else:
        left = ''; right = ''
    subgroups = repr_subgroups(subgroup(items, subgroup_size))
    dlim = ",\n" + " " * N_spaces
    return start + left + dlim.join(subgroups) + right

# %% Chemical management

def remove_undefined_chemicals(data: dict, chemicals: Chemicals):
    for i in tuple(data):
        if i not in chemicals: del data[i]

def default_chemical_dict(dct, chemicals, g, l, s, n=None):
    if n is None: n = l
    for i in chemicals:
        ID = i.ID
        if ID not in dct:
            locked_state = i.locked_state
            if locked_state == 'g': dct[ID] = g
            elif locked_state == 'l': dct[ID] = l
            elif locked_state == 's': dct[ID] = s
            elif locked_state is None: dct[ID] = n
            else: raise RuntimeError(f"unknown locked state '{locked_state}' of chemical '{i}'")
