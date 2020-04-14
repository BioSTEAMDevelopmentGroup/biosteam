# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:13:07 2019

@author: yoelr
"""
from biosteam import units
import doctest

__all__ = ('test_binary_distillation',
           'test_mixer',
           'test_splitter',
           'test_liquids_centrifuge',
           'test_fermentation',
           'test_heat_exchanger',
           'test_pump',
           'test_tank',
           'test_molecular_sieve',
           'test_flash',
           'test_mass_balance',
           'test_process_specification',
           'test_unit_operations')

def test_binary_distillation():
    doctest.testmod(units._binary_distillation)

def test_shortcut_column():
    doctest.testmod(units._shortcut_column)
    
def test_mixer():
    doctest.testmod(units._mixer)

def test_splitter():
    doctest.testmod(units._splitter)
    
def test_liquids_centrifuge():
    doctest.testmod(units._liquids_centrifuge)
    
def test_tank():
    doctest.testmod(units._tank)
    
def test_fermentation():
    doctest.testmod(units._fermentation)

def test_heat_exchanger():
    doctest.testmod(units._hx)

def test_pump():
    doctest.testmod(units._pump)

def test_molecular_sieve():
    doctest.testmod(units._molecular_sieve)

def test_flash():
    doctest.testmod(units._flash)
    
def test_mass_balance():
    doctest.testmod(units._balance)
    
def test_process_specification():
    doctest.testmod(units._process_specification)
    
def test_unit_operations():
    test_binary_distillation()
    test_shortcut_column()
    test_mixer()
    test_splitter()
    test_liquids_centrifuge()
    test_tank()
    test_fermentation()
    test_heat_exchanger()
    test_pump()
    test_molecular_sieve()
    test_flash()
    test_mass_balance()
    test_process_specification()
    
if __name__ == '__main__':
    test_unit_operations()