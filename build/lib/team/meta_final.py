# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:48:29 2019

@author: yoelr
"""

from .unit import metaUnit

class metaFinal(metaUnit):
    """Make Unit class final/unable to be inherited."""
    def __new__(mcl, name, bases, dct):
        # Make this class a final class
        for b in bases:
            if isinstance(b, mcl):
                raise TypeError(f"Cannot inherit from {b}. Instances of {mcl.__name__} cannot be inherited.")

        return super().__new__(mcl, name, bases, dct)