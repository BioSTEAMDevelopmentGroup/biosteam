# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 00:35:01 2019

@author: yoelr
"""
__all__ = ('NotImplementedMethod',)

class NotImplementedMethodType:
    __slots__ = ()
    
    @property
    def __name__(self): return "NotImplementedMethod"
    @property
    def __doc__(self): return None
    def __new__(self): return NotImplementedMethod
    def __call__(self): pass
    def __bool__(self): return False
    def __repr__(self): return "NotImplementedMethod"

NotImplementedMethod = object.__new__(NotImplementedMethodType)