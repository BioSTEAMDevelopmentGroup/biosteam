# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:24:01 2018

This module includes arbitrary classes and functions.

@author: Guest Group
"""
import biosteam as bst

__all__ = ('factor', 'checkbounds', 'approx2step', 'strtuple',
           'Counter', 'format_unit_line', 'format_unit_name')

class Counter:
    __slots__ = ('msg', 'N')
    def __init__(self, msg=None, N=0):
        self.msg = msg
        self.N = N
        
    def notify(self):
        print(f"{self.msg or 'counter'}: {self.N}")
        
    def restart(self, msg=None, N=0, notify=True):
        if notify: self.notify()
        self.msg = msg or self.msg
        self.N = N
        
    def count(self):
        self.N += 1
        
    def __repr__(self):
        return f"<Counter: msg={repr(self.msg)}, N={self.N}>"

# %% Number functions

def factor(base_units, new_units):
    if base_units == new_units: return 1
    else: return bst._Q(1, base_units).to(new_units).magnitude

def checkbounds(x, bounds):
    return bounds[0] < x < bounds[1]

def approx2step(val, x0, dx):
    """Approximate value, val, to closest increment/step, dx, starting from x0."""
    while True:
        if x0 > val: break
        x0 += dx
    return x0


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
        
def format_unit_line(line):
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
    