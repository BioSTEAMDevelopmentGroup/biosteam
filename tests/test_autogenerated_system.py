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
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, M1, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()
    
    
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
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, P3, M1, M2, S2, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()

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
    recycle_loop_sys.flatten()
    assert recycle_loop_sys.path == (P1, P2, P3, M1, M2, S2, S3, S1)
    recycle_loop_sys.empty_recycles()
    recycle_loop_sys.simulate()

if __name__ == '__main__':
    test_simple_recycle_loop()
    test_two_recycle_loops_with_complete_overlap()
    test_two_recycle_loops_with_partial_overlap()