# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for the heat exchanger network facility.
"""
import warnings
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose

def build_system(N_columns=1):
    """Doctest system of HeatExchangerNetwork; `N_columns > 1` adds more
    ShortcutColumns so auxiliary heat exchangers have duplicate IDs."""
    bst.settings.set_thermo(['Water', 'Methanol', 'Glycerol'], cache=True)
    bst.main_flowsheet.set_flowsheet('test_hxn')
    units = []
    feeds = []
    for i in range(N_columns):
        feed = bst.Stream(f'feed{i}', flow=(8000, 100 * (i + 1), 25))
        D = bst.ShortcutColumn(f'D{i}', ins=feed, LHK=('Methanol', 'Water'),
                               y_top=0.99, x_bot=0.01, k=2, is_divided=True)
        H1 = bst.HXutility(f'D{i}_H1', ins=D.outs[1], T=300)
        H2 = bst.HXutility(f'D{i}_H2', ins=D.outs[0], T=300)
        units.extend([D, H1, H2])
        feeds.append(feed)
    feed2 = bst.Stream('feed_flash', flow=(10000, 1000, 10))
    F1 = bst.Flash('F1', ins=feed2, V=0.9, P=101325)
    HXN = bst.HeatExchangerNetwork('HXN', T_min_app=5.)
    sys = bst.System.from_units('sys', units=[*units, F1, HXN])
    return sys, HXN, feeds[0]

def network_results(HXN):
    return dict(
        heat=HXN.actual_heat_util_load,
        cool=HXN.actual_cool_util_load,
        Q=np.array([hx.Q for hx in HXN.new_HXs]),
        installed=HXN.installed_costs['Heat exchangers'],
    )

def assert_same_results(a, b, rtol=1e-6):
    for key in a:
        assert_allclose(a[key], b[key], rtol=rtol, err_msg=key)

def simulate_cached(sys, HXN):
    HXN.cache_network = True
    HXN_sys = HXN.HXN_sys
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        sys.simulate()
    assert HXN.HXN_sys is HXN_sys, 'cached network was not used'

def test_cache_network_matches_fresh_synthesis():
    sys, HXN, feed = build_system()
    sys.simulate()
    fresh = network_results(HXN)
    assert HXN.actual_heat_util_load < 0.9 * HXN.original_heat_util_load
    simulate_cached(sys, HXN)
    assert_same_results(network_results(HXN), fresh)

def test_cache_network_perturbed_feed():
    sys, HXN, feed = build_system()
    sys.simulate()
    feed.F_mass *= 1.01
    simulate_cached(sys, HXN)
    cached = network_results(HXN)
    HXN.cache_network = False
    sys.simulate()
    assert_same_results(cached, network_results(HXN), rtol=1e-3)

def test_cache_network_duplicate_IDs():
    sys, HXN, feed = build_system(N_columns=2)
    sys.simulate()
    IDs = [hx.ID for hx in HXN.original_heat_exchangers]
    assert len(IDs) != len(set(IDs)), 'test needs duplicate auxiliary IDs'
    fresh = network_results(HXN)
    simulate_cached(sys, HXN)
    assert_same_results(network_results(HXN), fresh)

def test_energy_balance_error_contributions_ignored_none():
    sys, HXN, feed = build_system()
    sys.simulate()
    N = len(HXN.original_heat_utils)
    errors = HXN._energy_balance_error_contributions()
    assert len(errors) == N
    assert HXN.ignored is None

if __name__ == '__main__':
    test_cache_network_matches_fresh_synthesis()
    test_cache_network_perturbed_feed()
    test_cache_network_duplicate_IDs()
    test_energy_balance_error_contributions_ignored_none()
