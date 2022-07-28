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


def test_isentropic_hydrogen_compressor():
    bst.settings.set_thermo(["H2"])
    feed = bst.Stream(H2=1, T=25 + 273.15, P=101325, phase='g')
    K = bst.units.IsentropicCompressor(ins=feed, P=50e5, eta=0.7)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert allclose(
        a=[out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        b=[1.0, 0.0, 1151.3251608356125, 50e5],
    )
    # check compressor design
    assert allclose(
        a=list(K.design_results.values())[1:],
        b=[7.02839886238921, 1151.3251608356125, 24.465403697038127, 4.919879203671207, 901.1332666056242],
    )
    assert K.design_results["Type"] == "Blower"
    pass


def test_isentropic_two_phase_steam_compressor():
    bst.settings.set_thermo(["H2O"])
    # feed is steam-water-mixture
    feed = bst.MultiStream(T=372.75, P=1e5, l=[('H2O', 0.1)], g=[('H2O', 0.9)])
    # from compress until superheated steam (test case for vle parameter)
    K = bst.units.IsentropicCompressor(ins=feed, P=100e5, eta=1.0, vle=True)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert allclose(
        a=[out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        b=[1.0, 0.0, 797.7528062886108, 100e5],
    )
    # check compressor design
    assert allclose(
        a=list(K.design_results.values())[1:],
        b=[5.410389965766295, 797.7528062886108, 27.89482592365777, 5.410389965201038, 797.7528062360269],
    )
    assert K.design_results["Type"] == "Blower"
    pass


def test_isothermal_hydrogen_compressor():
    bst.settings.set_thermo(["H2"])
    feed = bst.Stream(H2=1, T=298.15, P=20e5, phase='g')
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert allclose(
        a=[out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        b=[1.0, 0.0, feed.T, 350e5],
    )
    # check compressor design
    assert allclose(
        a=list(K.design_results.values())[1:],
        b=[1.970908496944548, 298.15, 1.2394785148011942],
    )
    assert K.design_results["Type"] == "Blower"
    # check heat utility
    assert allclose(
        a=[K.heat_utilities[0].unit_duty, K.heat_utilities[0].duty, K.heat_utilities[0].flow],
        b=[-7095.270589000373, -7095.270589000373, 7.35376825845714],
    )
    pass


def test_compressor_design():
    bst.settings.set_thermo(["H2"])
    feed = bst.Stream(H2=1, T=298.15, P=20e5, phase='g')
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5)
    K.simulate()
    assert K.design_results["Type"] == "Blower"
    feed.F_mol = 100
    K.simulate()
    assert K.design_results["Type"] == "Reciprocating"
    feed.F_mol = 1e4
    K.simulate()
    assert K.design_results["Type"] == "Centrifugal"
    # test out of range
    with pytest.raises(RuntimeError):
        feed.F_mol = 2e4
        K.simulate()
    # test user overwrite
    feed.F_mol = 1e4
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5, type="Reciprocating")
    K.simulate()
    assert K.design_results["Type"] == "Reciprocating"
    pass


if __name__ == '__main__':
    test_compressor_design()
    test_isentropic_hydrogen_compressor()
    test_isentropic_two_phase_steam_compressor()
    test_isothermal_hydrogen_compressor()
