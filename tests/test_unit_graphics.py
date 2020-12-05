# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:25:23 2020

@author: yrc2
"""
import pytest
import biosteam as bst


def test_unit_graphics():
    sink = 'M1'
    sink_index = 2
    assert bst.Mixer._graphics.get_inlet_options(sink, sink_index) == {'headport': 'c'}
    source = 'M1'
    source_index = 0
    assert bst.Mixer._graphics.get_outlet_options(source, source_index) == {'tailport': 'e'}
    
    GraphicsWarning = bst.exceptions.GraphicsWarning
    with pytest.warns(GraphicsWarning):
        bst.Mixer._graphics.get_inlet_options(sink, 100)
    
    with pytest.warns(GraphicsWarning):
        bst.Mixer._graphics.get_outlet_options(sink, 1)
    