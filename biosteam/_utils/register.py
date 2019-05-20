# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

@author: Yoel Cortes-Pena
"""
from weakref import ref

__all__ = ('Register',)

_get = object.__getattribute__

def reload(self):
    d = dictionaries[self]
    for i, j in list(d.items()):
        j = j()
        if not j or i!=str(j): del d[i]
    must_reload.remove(self)
    
must_reload = set()
dictionaries = {}

class Register:
    
    def __init__(self):
        dictionaries[self] = _get(self, '__dict__')
    
    def __getattribute__(self, key):
        if self in must_reload: reload(self)
        attr = _get(self, key)
        if isinstance(attr, ref): return attr()
        return attr
    
    def __setattr__(self, key, value):
        dictionaries[self][key] = ref(value)
        must_reload.add(self)
    __setitem__ = __setattr__
    
    def __getitem__(self, key):
        if self in must_reload: reload(self)
        return _get(self, key)()
    
    def __iter__(self):
        if self in must_reload: reload(self)
        for i in dictionaries[self].values(): yield i()
    
    def __delitem__(self, key):
        dict = dictionaries[self]
        value = dict[key]
        if hasattr(value, '_delete'): value._delete()
        del dict[key]
    __delattr__ = __delitem__
    
    def __repr__(self):
        if self in must_reload: reload(self)
        return f'Register:\n ' + ',\n '.join([repr(i) for i in self])
    

    