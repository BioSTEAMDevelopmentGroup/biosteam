# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import numpy as np
from numpy.testing import assert_allclose
from biosteam import (
    main_flowsheet as f,
    Pump, Mixer, Splitter, HXprocess,
    Stream, settings,
    System, Network,
)

def test_simple_recycle_loop():
    f.set_flowsheet('simple_recycle_loop')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    product = Stream('product')
    P1 = Pump('P1', feedstock)
    P2 = Pump('P2', water)
    M1 = Mixer('M1', [P1-0, P2-0, recycle])
    S1 = Splitter('S1', M1-0, [product, recycle], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys.to_network()
    actual_network = Network(
        [P1,
         P2,
         Network(
            [M1,
             S1],
            recycle=recycle)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = recycle.mol.copy()
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, M1, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = recycle.mol.copy()
    assert_allclose(x_nested_solution, x_flat_solution, rtol=1e-2)
    
def test_two_recycle_loops_with_complete_overlap():
    f.set_flowsheet('two_recycle_loops_with_complete_overlap')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    inner_recycle = Stream('inner_recycle')
    product = Stream('product')
    inner_water = Stream('inner_water', Water=10)
    P1 = Pump('P1', feedstock)
    P2 = Pump('P2', water)
    P3 = Pump('P3', inner_water)
    M1 = Mixer('M1', [P1-0, P2-0, recycle])
    M2 = Mixer('M2', [M1-0, P3-0, inner_recycle])
    S2 = Splitter('S2', M2-0, ['', inner_recycle], split=0.5)
    S1 = Splitter('S1', S2-0, [product, recycle], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys.to_network()
    actual_network = Network(
        [P1,
         P2,
         P3,
         Network(
            [M1,
             Network(
                 [M2,
                  S2],
                 recycle=inner_recycle),
             S1],
            recycle=recycle)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle.mol, inner_recycle.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, P3, M1, M2, S2, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle.mol, inner_recycle.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=1e-2)

def test_two_recycle_loops_with_partial_overlap():
    f.set_flowsheet('two_recycle_loops_with_partial_overlap')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    inner_recycle = Stream('inner_recycle')
    product = Stream('product')
    coproduct = Stream('coproduct')
    inner_water = Stream('inner_water', Water=10)
    P1 = Pump('P1', feedstock)
    P2 = Pump('P2', water)
    P3 = Pump('P3', inner_water)
    M1 = Mixer('M1', [P1-0, P2-0, recycle])
    M2 = Mixer('M2', [M1-0, P3-0, inner_recycle])
    S2 = Splitter('S2', M2-0, split=0.5)
    S3 = Splitter('S3', S2-0, [coproduct, inner_recycle], split=0.5)
    S1 = Splitter('S1', S2-1, [product, recycle], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys.to_network()
    actual_network = Network(
        [P1,
         P2,
         P3,
         Network(
            [M1,
             Network(
                 [M2,
                  S2,
                  S3],
                 recycle=inner_recycle),
             S1],
            recycle=recycle)])
    actual_network_alternative = Network(
        [P1,
         P2,
         P3,
         Network(
            [Network(
                [M1,
                 M2,
                 S2,
                 S1],
                recycle=S1-1),
             S3],
            recycle=S3-1)])
    assert network == actual_network or network == actual_network_alternative
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle.mol, inner_recycle.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, P3, M1, M2, S2, S3, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle.mol, inner_recycle.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=1e-2)

def test_feed_forward_recycle_loop():
    f.set_flowsheet('feed_forward_recycle_loop')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    inner_recycle = Stream('inner_recycle')
    product = Stream('product')
    P1 = Pump('P1', feedstock)
    P2 = Pump('P2', water)
    M1 = Mixer('M1', [P1-0, recycle, inner_recycle])
    S1 = Splitter('S1', M1-0, ['', recycle], split=0.5)
    M2 = Mixer('M2', [P2-0, S1-0])
    S2 = Splitter('S2', M2-0, [inner_recycle, product], split=0.5)
    P3 = Pump('P3', product)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys.to_network()
    actual_network = Network(
        [P1,
         P2,
         Network(
            [M1,
             S1,
             M2,
             S2],
            recycle=S2-0),
         P3])
    assert network == actual_network 
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle.mol, inner_recycle.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, M1, S1, M2, S2, P3)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle.mol, inner_recycle.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=1e-2)

def test_nested_recycle_loops():
    f.set_flowsheet('feed_forward_recycle_loop')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    recycle_1, recycle_2, recycle_3, recycle_4, recycle_5 = recycles = [
        Stream('recycle_1'),
        Stream('recycle_2'),
        Stream('recycle_3'),
        Stream('recycle_4'),
        Stream('recycle_5'),
    ]
    feed_1 = Stream('feed_1', Water=10)
    feed_2 = Stream('feed_2', Water=10)
    feed_3 = Stream('feed_3', Water=10)
    feed_4 = Stream('feed_4', Water=10)
    inner_recycle = Stream('inner_recycle')
    product = Stream('product')
    P1 = Pump('P1', feedstock)
    M1 = Mixer('M1', [P1-0, recycle_1])
    P2 = Pump('P2', M1-0)
    H1 = HXprocess('H1', [P2-0, recycle_2])
    P3 = Pump('P3', feed_1)
    M3 = Mixer('M3', [H1-0, P3-0], recycle_2)
    P4 = Pump('P4', H1-1)
    M4 = Mixer('M4', [P4-0, recycle_3])
    P5 = Pump('P5', feed_2)
    P6 = Pump('P6', feed_3)
    M5 = Mixer('M5', [M4-0, P5-0, P6-0])
    H2 = HXprocess('H2', [M5-0, recycle_4])
    P7 = Pump('P7', feed_4)
    M6 = Mixer('M6', [H2-0, P7-0])
    S1 = Splitter('S1', M6-0, [recycle_4, ''], split=0.5)
    S2 = Splitter('S2', H2-1, split=0.5)
    P8 = Pump('P8', S2-1)
    S3 = Splitter('S3', P8-0, split=0.5)
    H3 = HXprocess('H3', [S2-0, recycle_5], ['', recycle_3])
    M7 = Mixer('M7', [H3-0, S3-0])
    S4 = Splitter('S4', M7-0, ['', recycle_5], split=0.5)
    S5 = Splitter('S5', S4-0, ['', recycle_1], split=0.5)
    M8 = Mixer('M8', [S3-1, S1-1, S5-0], product)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    recycle_loop_sys.simulate()
    network = recycle_loop_sys.to_network()
    actual_network = Network(
        [P1,
         P3,
         P5,
         P6,
         P7,
         Network(
            [M1,
             P2,
             Network(
                [H1,
                 M3],
                recycle=M3-0),
             P4,
             Network(
                [M4,
                 M5,
                 Network(
                    [H2,
                     M6,
                     S1],
                    recycle=S1-0),
                 S2,
                 Network(
                    [H3,
                     M7,
                     S4],
                    recycle=H3-0),
                 P8,
                 S3],
                recycle=H3-1),
             S5],
            recycle=S5-1),
         M8])
    assert network == actual_network 
    x_nested_solution = np.vstack([i.mol for i in recycles])
    recycle_loop_sys.flatten()
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([i.mol for i in recycles])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=1e-2)

if __name__ == '__main__':
    test_simple_recycle_loop()
    test_two_recycle_loops_with_complete_overlap()
    test_two_recycle_loops_with_partial_overlap()
    test_feed_forward_recycle_loop()
    test_nested_recycle_loops()