# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression tests for solids separation units.

SolidsCentrifuge purchase costs follow Seider et al. (2017), Table 16.32
(CE = 567), per centrifuge as a function of its solids loading S [ton/hr]:

* continuous scroll solid bowl: 68,040 S^0.50, S = 2-40 ton/hr
* continuous reciprocating pusher: 170,100 S^0.30, S = 1-20 ton/hr

Loadings above the upper limit are split evenly across ceil(S/S_max)
centrifuges in parallel.
"""
import biosteam as bst
from numpy.testing import assert_allclose

kg_per_ton = 1 / 0.0011023 # Short tons

def create_centrifuge(solids_loading, centrifuge_type):
    bst.settings.set_thermo([
        'Water',
        bst.Chemical('Solid', default=True, phase='s', MW=1, search_db=False),
    ], cache=True)
    feed = bst.Stream(
        Water=3 * solids_loading * kg_per_ton + 1000.,
        Solid=solids_loading * kg_per_ton,
        units='kg/hr',
    )
    centrifuge = bst.SolidsCentrifuge(
        ins=feed, split=dict(Solid=0.98), centrifuge_type=centrifuge_type,
    )
    centrifuge.simulate()
    return centrifuge

def test_solids_centrifuge_cost_depends_on_type():
    f = bst.CE / 567
    scroll = create_centrifuge(10., 'scroll_solid_bowl')
    pusher = create_centrifuge(10., 'reciprocating_pusher')
    assert scroll.design_results['Number of centrifuges'] == 1
    assert pusher.design_results['Number of centrifuges'] == 1
    assert_allclose(
        scroll.purchase_costs['Centrifuges'], 68040 * 10**0.5 * f, rtol=1e-6
    )
    assert_allclose(
        pusher.purchase_costs['Centrifuges'], 170100 * 10**0.3 * f, rtol=1e-6
    )

def test_solids_centrifuge_cost_scales_with_number_of_centrifuges():
    f = bst.CE / 567
    scroll = create_centrifuge(60., 'scroll_solid_bowl')
    pusher = create_centrifuge(60., 'reciprocating_pusher')
    assert scroll.design_results['Number of centrifuges'] == 2
    assert pusher.design_results['Number of centrifuges'] == 3
    scroll_cost = 2 * 68040 * 30**0.5 * f
    pusher_cost = 3 * 170100 * 20**0.3 * f
    assert_allclose(scroll.purchase_costs['Centrifuges'], scroll_cost, rtol=1e-6)
    assert_allclose(pusher.purchase_costs['Centrifuges'], pusher_cost, rtol=1e-6)
    assert_allclose(scroll.installed_costs['Centrifuges'], 2.03 * scroll_cost, rtol=1e-6)
    assert_allclose(pusher.installed_costs['Centrifuges'], 2.03 * pusher_cost, rtol=1e-6)

def test_solids_centrifuge_power_is_not_scaled_by_number_of_centrifuges():
    for centrifuge_type in ('scroll_solid_bowl', 'reciprocating_pusher'):
        centrifuge = create_centrifuge(60., centrifuge_type)
        assert centrifuge.design_results['Number of centrifuges'] > 1
        assert_allclose(
            centrifuge.power_utility.rate,
            centrifuge.kWhr_per_m3 * centrifuge.F_vol_in,
            rtol=1e-6,
        )

def test_solids_centrifuge_without_solids_has_no_cost():
    for centrifuge_type in ('scroll_solid_bowl', 'reciprocating_pusher'):
        centrifuge = create_centrifuge(0., centrifuge_type)
        assert centrifuge.design_results['Number of centrifuges'] == 0
        assert centrifuge.purchase_costs['Centrifuges'] == 0.

if __name__ == '__main__':
    test_solids_centrifuge_cost_depends_on_type()
    test_solids_centrifuge_cost_scales_with_number_of_centrifuges()
    test_solids_centrifuge_power_is_not_scaled_by_number_of_centrifuges()
    test_solids_centrifuge_without_solids_has_no_cost()
