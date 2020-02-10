# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:13:07 2019

@author: yoelr
"""
import biosteam as bst
import doctest

__all__ = ('test_binary_distillation', 'test_mixer', 'test_splitter',
           'test_tank')

def test_binary_distillation():
    doctest.testmod(bst.units._distillation)
    
def test_mixer():
    doctest.testmod(bst.units._mixer)

def test_splitter():
    doctest.testmod(bst.units._splitter)
    
def test_tank():
    doctest.testmod(bst.units._tank)
    
def test_unit_operations():
    test_binary_distillation()
    test_mixer()
    test_splitter()
    test_tank()
    
if __name__ == '__main__':
    test_unit_operations()