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
import numpy as np
from numpy.testing import assert_allclose

def test_compressor_drivers():
    bst.F.set_flowsheet('compressor_drivers')
    bst.settings.set_thermo(['N2', 'O2'], cache=True)
    feed = bst.Stream(N2=79, O2=21, T=300, phase='g')
    K = bst.IsentropicCompressor(ins=feed, P=5e5, eta=0.8, driver='Steam turbine')
    K.simulate()
    assert K.outs[0].T > feed.T
    assert K.power_utility.production > 0
    assert_allclose(K.power_utility.rate, 0, atol=1e-6)
    assert K.heat_utilities[0].duty > 0
    assert K.purchase_cost > 0
    K.driver = 'Gas turbine'
    with pytest.raises(RuntimeError):
        K.simulate()

def test_compressor_high_pressure_and_errors():
    bst.F.set_flowsheet('compressor_errors')
    bst.settings.set_thermo(['N2', 'O2'], cache=True)
    feed = bst.Stream(N2=79, O2=21, T=300, phase='g')
    K = bst.IsentropicCompressor(ins=feed, P=1e10, eta=0.8)
    with pytest.warns(bst.exceptions.DesignWarning):
        K.simulate()
    K = bst.IsothermalCompressor(ins=feed, P=5e4, eta=0.8)
    with pytest.raises(RuntimeError):
        K.simulate() # Outlet pressure lower than inlet
    with pytest.raises(ValueError):
        bst.PolytropicCompressor(ins=feed, P=5e5, method='not a method')

def test_isentropic_compressor_phenomena():
    bst.F.set_flowsheet('compressor_phenomena')
    bst.settings.set_thermo(['N2', 'O2'], cache=True)
    feed = bst.Stream(N2=79, O2=21, T=300, phase='g')
    K = bst.IsentropicCompressor(ins=feed, P=5e5, eta=0.8)
    K.simulate()
    assert K.T_node is K.E_node
    assert K.get_E_node(feed) is K.T_node
    assert K.Q_node is K.Q_node
    Q = K.Q
    K._update_nonlinearities()
    assert_allclose(K.Q, Q, rtol=1e-2)
    equations = K._create_energy_balance_equations()
    assert len(equations) == 1

def test_multistage_compressor_errors():
    bst.F.set_flowsheet('multistage_compressor_errors')
    bst.settings.set_thermo(['N2', 'O2'], cache=True)
    feed = bst.Stream(N2=79, O2=21, T=300, phase='g')
    ks = [bst.IsentropicCompressor(P=2e5), bst.IsentropicCompressor(P=4e5)]
    hxs = [bst.HXutility(T=300), bst.HXutility(T=300)]
    with pytest.raises(RuntimeError):
        bst.MultistageCompressor(ins=feed, compressors=[None], hxs=hxs)
    with pytest.raises(RuntimeError):
        bst.MultistageCompressor(ins=feed, compressors=ks, hxs=[None])
    with pytest.raises(RuntimeError):
        bst.MultistageCompressor(ins=feed, compressors=ks, hxs=hxs[:1])
    with pytest.raises(RuntimeError):
        bst.MultistageCompressor(ins=feed)
    K = bst.MultistageCompressor(ins=feed, pr=2, n_stages=2, eta=0.8)
    K.simulate()
    assert_allclose(K.outs[0].P, 4 * feed.P, rtol=1e-6)
    K.reset_cache()
    K.simulate()
    assert K.purchase_cost > 0

if __name__ == '__main__':
    test_compressor_drivers()
    test_compressor_high_pressure_and_errors()
    test_isentropic_compressor_phenomena()
    test_multistage_compressor_errors()
