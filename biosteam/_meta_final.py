# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:48:29 2019

@author: yoelr
"""

from ._unit import metaUnit

class metaFinal(metaUnit):
    """Make Unit class final/unable to be inherited."""
    def __new__(mcl, name, bases, dct):
        try:
            base, = bases
        except:
            TypeError('cannot create {mcl.__name__} Unit subclass from more than one super class')
            
        if isinstance(base, mcl):
            raise TypeError(f"cannot inherit from {base}. Instances of {mcl.__name__} cannot be inherited")

        return super().__new__(mcl, name, bases, dct)