# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 16:47:33 2018

@author: Yoel Cortes-Pena
"""
from weakref import ref

__all__ = ('Register', 'search_register')

getattr_ = object.__getattribute__
setattr_ = object.__setattr__

def _reload(self):
    dict_ = getattr_(self, '__dict__')
    # Dictionary will change size in for loop, so use tuple
    for i, r in tuple(dict_.items()):
        j = r()
        if not j:
            del dict_[i]
        else:
            s = str(j)
            if i!=s: del dict_[i]
            dict_[s] = r
    _must_reload.remove(self)

do_nothing = lambda: None

def search_register(self, key):
    if self in _must_reload: _reload(self)
    return getattr_(self, '__dict__').get(key, do_nothing)()

_must_reload = set()

class Register:
    
    def __getattribute__(self, key):
        if self in _must_reload: _reload(self)
        try: return getattr_(self, key)()
        except: return getattr_(self, key)
    
    __getitem__ = __getattribute__
    
    def __bool__(self):
        return bool(getattr_(self, '__dict__'))
    
    def __setattr__(self, key, value):
        """Register object."""
        ID_words = key.split('_')
        assert all(word.isalnum() for word in ID_words), (
                'ID may only contain letters, numbers, and/or underscores; '
                'no special characters or spaces')
        value._ID = key
        setattr_(self, key, ref(value))
        _must_reload.add(self)
   
    def __setitem__(self, key, value):
        setattr_(self, key, ref(value))
        _must_reload.add(self)
    
    def __iter__(self):
        if self in _must_reload: _reload(self)
        for i in getattr_(self, '__dict__').values(): yield i()
    
    def __delattr__(self, key):
        dict_ = getattr_(self, '__dict__')
        value = dict_[key]()
        if hasattr(value, '_disconnect'): value._disconnect()
        del dict_[key]
    __delitem__ = __delattr__
    
    def __repr__(self):
        if getattr_(self, '__dict__'):
            if self in _must_reload: _reload(self)
            return f'Register:\n ' + '\n '.join([repr(i) for i in self])
        else:
            return f'Register: (Empty)'

