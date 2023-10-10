# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import numpy as np
import biosteam as bst
from biosteam._network import Network
from numpy.testing import assert_allclose
from biosteam import (
    main_flowsheet as f,
    Pump, Mixer, Splitter, HXprocess,
    Stream, settings, System
)

def test_trivial_case():
    f.set_flowsheet('trivial_case')
    settings.set_thermo(['Water'], cache=True)
    trivial_sys_a = f.create_system('trivial_sys_a')
    network_a = trivial_sys_a._to_network()
    Stream('feedstock', Water=1000)
    Stream('water', Water=10)
    trivial_sys_b = f.create_system('trivial_sys_b')
    network_b = trivial_sys_b._to_network()
    actual_network = Network(
        [])
    assert network_a == network_b == actual_network
    trivial_sys_b.simulate()
    trivial_sys_b.flatten()
    assert trivial_sys_b.path == ()
    trivial_sys_b.empty_recycles()
    trivial_sys_b.simulate()
    f.clear()

def test_linear_case():
    f.set_flowsheet('linear_case')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    byproduct = Stream('byproduct')
    product = Stream('product')
    M1 = Mixer('M1', [feedstock, water])
    S1 = Splitter('S1', M1-0, [product, byproduct], split=0.5)
    linear_sys = f.create_system('linear_sys')
    network = linear_sys._to_network()
    actual_network = Network(
        [M1,
         S1])
    assert network == actual_network
    linear_sys.simulate()
    linear_sys.flatten()
    assert linear_sys.path == (M1, S1)
    linear_sys.empty_recycles()
    linear_sys.simulate()
    f.clear()

def test_simple_recycle_loop():
    f.set_flowsheet('simple_recycle_loop')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    product = Stream('product')
    M1 = Mixer('M1', [feedstock, water, recycle])
    S1 = Splitter('S1', M1-0, [product, recycle], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [M1,
         S1],
        recycle=recycle)
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = recycle.mol.copy()
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (M1, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = recycle.mol.copy()
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

def test_unconnected_case():
    f.set_flowsheet('unconnected_case')
    settings.set_thermo(['Water'], cache=True)
    # Recycle loop a
    feedstock_a = Stream('feedstock_a', Water=1000)
    water_a = Stream('water_a', Water=10)
    byproduct_a = Stream('byproduct_a')
    product_a = Stream('product_a')
    M1_a = Mixer('M1_a', [feedstock_a, water_a])
    S1_a = Splitter('S1_a', M1_a-0, [product_a, byproduct_a], split=0.5)
    # Recycle loop b
    feedstock_b = Stream('feedstock_b', Water=1000)
    water_b = Stream('water_b', Water=10)
    byproduct_b = Stream('byproduct_b')
    product_b = Stream('product_b')
    M1_b = Mixer('M1_b', [feedstock_b, water_b])
    S1_b = Splitter('S1_b', M1_b-0, [product_b, byproduct_b], split=0.5)
    parallel_sys = f.create_system('parallel_sys')
    network = parallel_sys._to_network()
    actual_network = Network(
        [M1_a,
         S1_a,
         M1_b,
         S1_b])
    assert network == actual_network
    parallel_sys.simulate()
    parallel_sys.flatten()
    assert parallel_sys.path == (M1_a, S1_a, M1_b, S1_b)
    parallel_sys.empty_recycles()
    parallel_sys.simulate()
    # Process specification reordering
    assert parallel_sys.path == (M1_a, S1_a, M1_b, S1_b)
    M1_b.add_specification(f=lambda: None, run=True, impacted_units=[M1_a])
    parallel_sys.simulate() # Should reorder due to process specifications
    assert parallel_sys.path == (M1_b, M1_a, S1_a, S1_b)
    f.clear()

def test_unconnected_recycle_loop():
    f.set_flowsheet('unconnected_recycle_loop')
    settings.set_thermo(['Water'], cache=True)
    # Recycle loop a
    feedstock_a = Stream('feedstock_a', Water=1000)
    water_a = Stream('water_a', Water=10)
    recycle_a = Stream('recycle_a')
    product_a = Stream('product_a')
    M1_a = Mixer('M1_a', [feedstock_a, water_a, recycle_a])
    S1_a = Splitter('S1_a', M1_a-0, [product_a, recycle_a], split=0.5)
    # Recycle loop b
    feedstock_b = Stream('feedstock_b', Water=800)
    water_b = Stream('water_b', Water=10)
    byproduct_b = Stream('byproduct_b')
    product_b = Stream('product_b')
    M1_b = Mixer('M1_b', [feedstock_b, water_b])
    S1_b = Splitter('S1_b', M1_b-0, [product_b, byproduct_b], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [Network(
            [M1_a,
             S1_a],
            recycle=S1_a-1),
         M1_b,
         S1_b])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle_a.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (M1_a, S1_a, M1_b, S1_b)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle_a.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    # New configuration due to stream priority changes
    feedstock_b.F_mol = 1200 # Larger flow means b has higher priority
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [M1_b,
         S1_b,
         Network(
            [M1_a,
             S1_a],
            recycle=S1_a-1)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle_a.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (M1_b, S1_b, M1_a, S1_a)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle_a.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    # Process specification reordering
    assert recycle_loop_sys.path == (M1_b, S1_b, M1_a, S1_a)
    M1_a.add_specification(f=lambda: None, run=True, impacted_units=[M1_b])
    recycle_loop_sys.simulate() # Should reorder due to process specifications
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (M1_a, S1_a, M1_b, S1_b)
    f.clear()

def test_unconnected_recycle_loops():
    f.set_flowsheet('unconnected_recycle_loops')
    settings.set_thermo(['Water'], cache=True)
    feedstock_a = Stream('feedstock_a', Water=1000)
    water_a = Stream('water_a', Water=10)
    recycle_a = Stream('recycle_a')
    product_a = Stream('product_a')
    M1_a = Mixer('M1_a', [feedstock_a, water_a, recycle_a])
    S1_a = Splitter('S1_a', M1_a-0, [product_a, recycle_a], split=0.5)
    
    feedstock_b = Stream('feedstock_b', Water=800)
    water_b = Stream('water_b', Water=10)
    recycle_b = Stream('recycle_b')
    product_b = Stream('product_b')
    M1_b = Mixer('M1_b', [feedstock_b, water_b, recycle_b])
    S1_b = Splitter('S1_b', M1_b-0, [product_b, recycle_b], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [Network(
            [M1_a,
             S1_a],
            recycle=S1_a-1),
         Network(
            [M1_b,
             S1_b],
            recycle=S1_b-1)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (M1_a, S1_a, M1_b, S1_b)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()
    
def test_inner_recycle_loop():
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
    network = recycle_loop_sys._to_network()
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
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()
    
def test_inner_recycle_loop_with_bifurcated_feed():
    f.set_flowsheet('simple_recycle_loop_with_bifurcated_feed')
    settings.set_thermo(['Water'], cache=True)
    feedstock = Stream('feedstock', Water=1000)
    water = Stream('water', Water=10)
    recycle = Stream('recycle')
    product = Stream('product')
    side_product = Stream('side_product')
    P1 = Pump('P1', feedstock)
    P2 = Pump('P2', water)
    S1 = Splitter('S1', P2-0, split=0.5)
    M1 = Mixer('M1', [P1-0, S1-1, recycle])
    S2 = Splitter('S2', M1-0, [side_product, ''], split=0.5)
    M2 = Mixer('M2', [S1-0, S2-1])
    S3 = Splitter('S3', M2-0, [product, recycle], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [P1,
         P2,
         S1,
         Network(
            [M1,
             S2,
             M2,
             S3],
            recycle=S3-1)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = recycle.mol.copy()
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, S1, M1, S2, M2, S3)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = recycle.mol.copy()
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

def test_bifurcated_recycle_loops():
    f.set_flowsheet('bifurcated_recycle_loops')
    settings.set_thermo(['Water'], cache=True)
    feed_a = Stream('feed_a', Water=10)
    water_a = Stream('water_a', Water=10)
    recycle_a = Stream('recycle_a')
    P1_a = Pump('P1_a', feed_a)
    P2_a = Pump('P2_a', water_a)
    S1_a = Splitter('S1_a', P2_a-0, split=0.5)
    M1_a = Mixer('M1_a', [P1_a-0, S1_a-1, recycle_a])
    S2_a = Splitter('S2_a', M1_a-0, split=0.5)
    M2_a = Mixer('M2_a', [S1_a-0, S2_a-1])
    S3_a = Splitter('S3_a', M2_a-0, ['', recycle_a], split=0.5)
    
    feedstock_b = Stream('feedstock_b', Water=1000)
    recycle_b = Stream('recycle_b')
    product_b = Stream('product_b')
    side_product_b = Stream('side_product_b')
    P1_b = Pump('P1_b', feedstock_b)
    P2_b = Pump('P2_b', S2_a-0)
    S1_b = Splitter('S1_b', P2_b-0, split=0.5)
    M1_b = Mixer('M1_b', [P1_b-0, S1_b-1, recycle_b])
    S2_b = Splitter('S2_b', M1_b-0, [side_product_b, ''], split=0.5)
    M2_b = Mixer('M2_b', [S1_b-0, S2_b-1, S3_a-0])
    S3_b = Splitter('S3_b', M2_b-0, [product_b, recycle_b], split=0.5)
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    recycle_loop_sys.set_tolerance(mol=1e-3, rmol=1e-6, subsystems=True)
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [P1_b,
         P1_a,
         P2_a,
         S1_a,
         Network(
            [M1_a,
             S2_a,
             M2_a,
             S3_a],
            recycle=S3_a-1),
         P2_b,
         S1_b,
         Network(
            [M1_b,
             S2_b,
             M2_b,
             S3_b],
            recycle=S3_b-1)])
    assert network == actual_network
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    # Test flattend solution
    recycle_loop_sys.flatten()
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle_a.mol, recycle_b.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

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
    network = recycle_loop_sys._to_network()
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
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

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
    network = recycle_loop_sys._to_network()
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
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

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
    network = recycle_loop_sys._to_network()
    actual_network = (
        Network(
            [P1,
             P2,
             Network(
                [M1,
                 S1,
                 M2,
                 S2],
                recycle=M1-0),
             P3])
    )
    assert network == actual_network 
    recycle_loop_sys.simulate()
    x_nested_solution = np.vstack([recycle.mol, inner_recycle.mol])
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, M1, S1, M2, S2, P3)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([recycle.mol, inner_recycle.mol])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    f.clear()

def test_separate_recycle_loops():
    f.set_flowsheet('separate_recycle_loops')
    settings.set_thermo(['Water'], cache=True)
    feedstock_a = Stream('feedstock_a', Water=1000)
    water_a = Stream('water_a', Water=10)
    recycle_a = Stream('recycle_a')
    product_a = Stream('product_a')
    P1_a = Pump('P1_a', feedstock_a)
    P2_a = Pump('P2_a', water_a)
    M1_a = Mixer('M1_a', [P1_a-0, P2_a-0, recycle_a])
    S1_a = Splitter('S1_a', M1_a-0, [product_a, recycle_a], split=0.5)
    feedstock_b = Stream('feedstock_b', Water=1000)
    water_b = Stream('water_b', Water=10)
    recycle_b = Stream('recycle_b')
    product_b = Stream('product_b')
    P1_b = Pump('P1_b', feedstock_b)
    P2_b = Pump('P2_b', water_b)
    M1_b = Mixer('M1_b', [P1_b-0, P2_b-0, recycle_b])
    S1_b = Splitter('S1_b', M1_b-0, [product_b, recycle_b], split=0.5)
    recycles = [recycle_a, recycle_b]
    recycle_loop_sys = f.create_system('recycle_loop_sys')
    recycle_loop_sys.simulate()
    network = recycle_loop_sys._to_network()
    actual_network = Network(
        [P1_a,
         P2_a,
         Network(
            [M1_a,
             S1_a],
            recycle=S1_a-1),
         P1_b,
         P2_b,
         Network(
            [M1_b,
             S1_b],
            recycle=S1_b-1)])
    assert network == actual_network
    x_nested_solution = np.vstack([i.mol for i in recycles])
    recycle_loop_sys.flatten()
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    x_flat_solution = np.vstack([i.mol for i in recycles])
    assert_allclose(x_nested_solution, x_flat_solution, rtol=2e-2)
    assert recycle_loop_sys.path == (P1_a, P2_a, M1_a, S1_a, P1_b, P2_b, M1_b, S1_b)
    f.clear()
    
def test_nested_recycle_loops():
    f.set_flowsheet('nested_recycle_loops')
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
    network = recycle_loop_sys._to_network()
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
                 P8,
                 S3,
                 Network(
                    [M7,
                     S4,
                     H3],
                    recycle=H3-0)],
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
    assert_allclose(x_nested_solution, x_flat_solution, rtol=5e-2)
    f.clear()

def test_sugarcane_ethanol_biorefinery_network():
    from biorefineries import sugarcane as sc
    sugarcane_sys = bst.System.from_units('sugarcane_sys', sc.sys.units)
    network = sugarcane_sys._to_network()
    actual_network = Network(
        [sc.U101,
         sc.U102,
         sc.U103,
         Network(
            [sc.U201,
             sc.S201,
             sc.M201],
            recycle=sc.M201-0),
         sc.T202,
         sc.H201,
         sc.T203,
         sc.P201,
         sc.T204,
         sc.T205,
         sc.P202,
         Network(
            [sc.M202,
             sc.H202,
             sc.T206,
             sc.C201,
             sc.C202,
             sc.P203],
            recycle=sc.P203-0),
         sc.S202,
         sc.F301,
         sc.P306,
         sc.M301,
         sc.H301,
         Network(
            [sc.R301,
             sc.T301,
             sc.C301,
             sc.S302],
            recycle=sc.S302-0),
         sc.D301,
         sc.M302,
         sc.P301,
         Network(
            [sc.H302,
             sc.D302,
             sc.P302],
            recycle=sc.P302-0),
         Network(
            [sc.M303,
             sc.D303,
             sc.H303,
             sc.U301],
            recycle=sc.U301-0),
         sc.H304,
         sc.T302,
         sc.P304,
         sc.T303,
         sc.P305,
         sc.M304,
         sc.T304,
         sc.P303,
         sc.M305,
         sc.U202,
         sc.M401])
    assert network == actual_network
    bst.process_tools.default()
    f.clear()
    
if __name__ == '__main__':
    test_trivial_case()
    test_linear_case()
    test_unconnected_case()
    test_simple_recycle_loop()
    test_unconnected_recycle_loop()
    test_unconnected_recycle_loops()
    test_inner_recycle_loop()
    test_inner_recycle_loop_with_bifurcated_feed()
    test_bifurcated_recycle_loops()
    test_two_recycle_loops_with_complete_overlap()
    test_two_recycle_loops_with_partial_overlap()
    test_feed_forward_recycle_loop()
    test_separate_recycle_loops()
    test_nested_recycle_loops()
    test_sugarcane_ethanol_biorefinery_network()