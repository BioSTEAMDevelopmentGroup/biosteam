# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:25:23 2020

@author: yrc2
"""
import pytest
import biosteam as bst

def test_unit_connections():
    from biorefineries.sugarcane import flowsheet as f
    globals().update(f.unit.__dict__)
    assert R301.neighborhood(1) == {D301, H301, T301, T305}
    assert R301.neighborhood(2) == {C301, D301, H301, M301, M302, R301, T301, T305}
    assert R301.neighborhood(100) == R301.neighborhood(1000) == {P202, T201, U201, D301, M304, 
                                                                 P201, P303, S201, T304, M201, 
                                                                 H303, S301, C201, P301, H301, 
                                                                 C202, U301, T205, P302, H302, 
                                                                 T305, T301, M303, C301, H304, 
                                                                 T204, M202, S202, P306, P203, 
                                                                 BT,   U202, T302, U101, F301, 
                                                                 D303, U102, M302, H202, H201, 
                                                                 P304, U103, T202, M305, T206, 
                                                                 M301, D302, T303, R301, T203, 
                                                                 P305}
    assert R301.get_downstream_units() == {T301, M303, H303, P304, C301, 
                                           H304, P301, M305, D301, T302, 
                                           M304, D302, U301, D303, P302, 
                                           H302, P303, M302, T304}
    
    ins = tuple(R301.ins)
    outs = tuple(R301.outs)
    R301.disconnect()
    assert not any(R301.ins + R301.outs)
    ins - R301 - outs
    
    class DummyUnit(bst.Unit, isabstract=True, new_graphics=False):
        _N_ins = 2; _N_outs = 2
        
    unit = DummyUnit(None)
    unit.take_place_of(R301)
    assert tuple(unit.ins + unit.outs) == (ins + outs)
    assert not any(R301.ins + R301.outs)
    R301.take_place_of(unit)
    
    
if __name__ == '__main__':
    test_unit_connections()