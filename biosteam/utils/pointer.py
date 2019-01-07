# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 18:57:13 2018

@author: Guest Group
"""

# Do not include '__new__', '__init__', '__repr__'
_magic_names = ('__del__', '__bytes__',
                '__str__', '__format__', '__lt__', 
                '__le__', '__eq__', '__ne__', '__gt__', '__ge__',
                '__hash__', '__bool__', '__getattr__', '__getattribute__',
                '__setattr__', '__delattr__', '__dir__', '__get__', '__set__',
                '__delete__', '__set_name__', '__call__', '__len__',
                '__getitem__', '__setitem__', '__delitem__', '__missing__',
                '__reversed__', '__contains__', '__add__', '__sub__', '__mul__',
                '__matmul__', '__truediv__', '__floordiv__', '__mod__',
                '__divmod__', '__pow__', '__lshift__', '__rshift__',
                '__and__', '__or__', '__radd__', '__rsub__', '__rmul__',
                '__rmatmul__', '__rtruediv__', '__rfloordiv__', '__rmod__',
                '__rdivmod__', '__rpow__', '__rlshift__', '__rrshift__',
                '__rand__', '__rxor__', '__ror__', '__iadd__', '__isub__',
                '__imul__', '__imatmul__', '__itruediv__', '__ifloordiv__',
                '__imod__', '__ipow__', '__ilshift__', '__irshift__', '__iand__',
                '__ixor__', '__ior__', '__neg__', '__pos__', '__abs__',
                '__invert__', '__complex__', '__int__', '__index__',
                '__round__', '__trunc__', '__floor__', '__ceil__',
                '__enter__', '__exit__', '__await__', '__aiter__', 
                '__anext__', '__aenter__', '__aexit__')


class metaPointer(type):
    
    @staticmethod
    def __wrap(name):
        def proxy_method(self, *args, **kwargs):
            value = type(self).get_value(self)
            method = getattr(value, name)
            return method(*args, **kwargs)
        return proxy_method
    
    def __new__(mcl, clsname, superclasses, new_definitions):
        wrap = mcl.__wrap
        for method_name in _magic_names:
            method = new_definitions.get(method_name)
            if not method:
                new_definitions[method_name] = wrap(method_name)
        return type.__new__(mcl, clsname, superclasses, new_definitions)


class Pointer:
    __slots__ = ('obj', 'key')
    
    def get_source(self):
        obj = Pointer.__getattribute__(self, 'obj')
        key = Pointer.__getattribute__(self, 'key')
        return obj, key
    
    def __init__(self, obj, key):
        super().__setattr__('obj', obj)
        super().__setattr__('key', key)
        
    def __repr__(self):
        value = type(self).get_value(self)
        return  f'{type(self).__name__}:{value}'


class Attribute(Pointer, metaclass=metaPointer):
    __slots__ = ('obj', 'key')
    
    def get_value(self):
        obj, key = Pointer.get_source(self)
        return getattr(obj, key)
    

class Item(Pointer, metaclass=metaPointer):
    __slots__ = ('obj', 'key')
    
    def get_value(self):
        obj, key = Pointer.get_source(self)
        return obj[key]
    