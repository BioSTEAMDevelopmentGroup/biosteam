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

__all__ = ('factor', 'checkbounds', 'strtuple',
           'format_title', 'format_unit_name',
           'remove_undefined_chemicals')

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
    for i in line:
        if i.isupper():
            words.append(word)
            word = i
        else:
            word += i
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
        if not all([i.isupper() for i in word]):
            word = word.lower()
        words.append(word)
    return ' '.join(words)

def format_unit_name(name):
    words = name.split(' ')
    new_words = []
    for i in words:
        new_words.append(i[0].capitalize() + i[1:])
    return ''.join(new_words)
    

# %% Chemical management

def remove_undefined_chemicals(data: dict, chemicals: Chemicals):
    for i in tuple(data):
        if i not in chemicals: del data[i]
