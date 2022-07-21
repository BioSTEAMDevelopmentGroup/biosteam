# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
from numpy import allclose

def test_hydrogen_compressor():
    bst.settings.set_thermo(["H2"])
    feed = bst.Stream(H2=1, T=25 + 273.15, P=101325, phase='g')
    C = bst.units.IsentropicCompressor(ins=feed, P=50e5, eta=0.7)
    C.simulate()
    assert allclose(
        a=list(C.design_results.values()),
        b=[7.02839886238921, 4.919879203671207, 1151.3251608356125, 901.1332666056242],
    )
    pass

    
if __name__ == '__main__':
    test_hydrogen_compressor()