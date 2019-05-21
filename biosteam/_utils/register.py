# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

@author: Yoel Cortes-Pena
"""
from weakref import ref

__all__ = ('Register', 'search_register')

_attr = object.__getattribute__

def _reload(self):
    dict_ = _dictionaries[self]
    for i, j in list(dict_.items()):
        j = j()
        if not j or i!=str(j): del dict_[i]
    _must_reload.remove(self)
    
_must_reload = set()
_dictionaries = {}

def search_register(self, key):
    ref = _dictionaries[self].get(key)
    return ref() if ref else ref 

class Register:
    
    def __init__(self):
        _dictionaries[self] = _attr(self, '__dict__')
    
    def __getattribute__(self, key):
        if self in _must_reload: _reload(self)
        attr = _attr(self, key)
        if isinstance(attr, ref): return attr()
        return attr
    
    def __setattr__(self, key, value):
        _dictionaries[self][key] = ref(value)
        _must_reload.add(self)
    __setitem__ = __setattr__
    
    def __getitem__(self, key):
        if self in _must_reload: _reload(self)
        return _attr(self, key)()
    
    def __iter__(self):
        if self in _must_reload: _reload(self)
        for i in _dictionaries[self].values(): yield i()
    
    def __delitem__(self, key):
        dict_ = _dictionaries[self]
        value = dict_[key]
        if hasattr(value, '_delete'): value._delete()
        del dict_[key]
    __delattr__ = __delitem__
    
    def __repr__(self):
        if self in _must_reload: _reload(self)
        return f'Register:\n ' + ',\n '.join([repr(i) for i in self])
    

    