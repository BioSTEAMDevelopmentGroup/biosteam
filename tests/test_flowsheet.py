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

def test_flowsheet_magic_methods():
    F = bst.F
    F.clear()
    flowsheet = F.get_flowsheet()
    assert flowsheet in F.flowsheet
    assert flowsheet.ID in F.flowsheet
    assert 12 not in F.flowsheet
    with pytest.raises(TypeError): 
        F.flowsheet.default_flowsheet = flowsheet
    assert len(list(F.flowsheet)) >= 1
    F.set_flowsheet('default')
    with pytest.raises(AttributeError):
        del F.flowsheet.default
    
    F.set_flowsheet('new_flowsheet')
    del F.flowsheet.default
    assert 'default' not in F.flowsheet
    
    with pytest.raises(AttributeError):
        F.new_flowsheet = F.new_flowsheet
    
def test_flowsheet_search():
    bst.settings.set_thermo(['Water'], cache=True)
    mixer = bst.Mixer()
    assert mixer in bst.F.unit
    bst.F.discard(mixer)
    assert mixer not in bst.F.unit
    
if __name__ == '__main__':
    test_flowsheet_magic_methods()
    test_flowsheet_search()