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
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.digraph import digraph

def create_system(ID):
    bst.F.set_flowsheet(ID)
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=1)
    other_feed = bst.Stream('other_feed', Water=5)
    empty_feed = bst.Stream('empty_feed')
    recycle = bst.Stream('recycle')
    M0 = bst.Mixer('M0', feed, 's0')
    M1 = bst.Mixer('M1', [M0-0, recycle, other_feed, empty_feed])
    S1 = bst.Splitter('S1', M1-0, ['s1', recycle], split=0.5)
    S2 = bst.Splitter('S2', S1-0, ['product', 'other_product'], split=0.5)
    S3 = bst.Splitter('S3', S2-1, ['empty_product', 'last_product'], split=0.)
    subsystem = bst.System('subsystem', path=[M1, S1], recycle=recycle)
    system = bst.System(ID, path=[M0, subsystem, S2, S3])
    system.simulate()
    return system

def test_system_diagram_kinds():
    system = create_system('system_diagram_kinds')
    preferences = bst.preferences
    with preferences.temporary():
        preferences.stream_width = 'F_mass'
        for kind in ('cluster', 'thorough', 'surface', 'minimal'):
            f = system.diagram(kind, display=False, number=True, label=True)
            assert 'digraph' in f.source
        preferences.show_all_streams = False
        preferences.fill_cluster = True
        for kind in ('cluster', 'surface', 'minimal'):
            f = system.diagram(kind, display=False)
            assert 'digraph' in f.source
    with pytest.raises(ValueError):
        system.diagram('invalid', display=False)

def test_surface_digraph_with_multiple_feeds_and_products():
    system = create_system('surface_digraph_with_multiple_feeds_and_products')
    F = system.flowsheet
    subsystem = bst.System('subsystem_2', path=[F.M0, F.M1, F.S1, F.S2, F.S3],
                           recycle=F.recycle)
    for show_all_streams in (True, False):
        with bst.preferences.temporary():
            bst.preferences.show_all_streams = show_all_streams
            f = digraph.surface_digraph([subsystem, subsystem])
            assert 'digraph' in f.source
    empty_1 = bst.Stream('empty_1')
    empty_2 = bst.Stream('empty_2')
    M2 = bst.Mixer('M2', [empty_1, empty_2], ['empty_3'])
    S4 = bst.Splitter('S4', M2-0, ['empty_4', 'empty_5'], split=0.5)
    subsystem = bst.System('subsystem_3', path=[M2, S4])
    with bst.preferences.temporary():
        bst.preferences.show_all_streams = False
        f = digraph.surface_digraph([subsystem])
        assert 'digraph' in f.source

def test_pen_width():
    bst.F.set_flowsheet('pen_width')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    streams = [bst.Stream(None, Water=i) for i in (1, 2, 3)]
    pen_width = digraph.PenWidth('F_mass', streams)
    widths = [float(pen_width(i)) for i in streams]
    assert widths == sorted(widths)
    big = bst.Stream(None, Water=100)
    assert float(pen_width(big)) >= widths[-1]

def test_fix_valve_symbol_in_svg_output():
    svg = (
        b'<svg xmlns="http://www.w3.org/2000/svg"><g><text x="0">V1</text>'
        b'<g><g><polygon points="0,0 1,0 1,1 0,1"/>'
        b'<image width="10px" x="5"/></g></g></g></svg>'
    )
    img = digraph.fix_valve_symbol_in_svg_output(
        svg, unit_color='white:#CDCDCD', unit_periphery_color='black',
        label_color='black',
    )
    assert b'image' not in img
    assert b'#CDCDCD' in img
    img = digraph.fix_valve_symbol_in_svg_output(svg)
    assert b'image' not in img

if __name__ == '__main__':
    test_system_diagram_kinds()
    test_surface_digraph_with_multiple_feeds_and_products()
    test_pen_width()
    test_fix_valve_symbol_in_svg_output()
