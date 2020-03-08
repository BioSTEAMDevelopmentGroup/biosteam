# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 06:45:34 2020

@author: yoelr
"""
import biosteam as bst
import doctest

def test_heat_utility():
    from biosteam import _heat_utility
    doctest.testmod(_heat_utility)
    
def test_power_utility():
    from biosteam import _power_utility
    doctest.testmod(_power_utility)
    
def test_utilities():
    test_heat_utility()
    test_power_utility()
    
if __name__ == '__main__':
    test_utilities()