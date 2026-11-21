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

def check_costs(unit):
    assert unit.purchase_cost > 0
    assert np.isfinite(unit.purchase_cost)
    assert unit.installed_cost >= unit.purchase_cost

def test_hxutility_exchanger_types():
    bst.F.set_flowsheet('hx_types')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    for kind in ('Fixed head', 'U tube', 'Kettle vaporizer', 'Floating head'):
        feed = bst.Stream(Water=2000, Ethanol=200, T=300)
        hx = bst.HXutility(ins=feed, T=350, rigorous=False, heat_exchanger_type=kind)
        hx.simulate()
        assert hx.heat_exchanger_type == kind
        assert kind in hx.baseline_purchase_costs
        check_costs(hx)
        assert hx.heat_utilities[0].duty > 0
    # Small area -> double pipe with area below minimum
    feed = bst.Stream(Water=0.01, T=300)
    hx = bst.HXutility(ins=feed, T=310)
    hx.simulate()
    assert 'Double pipe' in hx.baseline_purchase_costs

def test_hxutility_rigorous_specifications():
    bst.F.set_flowsheet('hx_rigorous')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, Ethanol=10, T=300)
    hx = bst.HXutility(ins=feed, V=0.5, rigorous=True)
    hx.simulate()
    assert_allclose(hx.outs[0].vapor_fraction, 0.5, rtol=1e-2)
    hx.V = 1
    hx.simulate()
    assert hx.outs[0].phase == 'g'
    assert hx.outs[0].T > feed.T
    hx.V = 0
    hx.simulate()
    assert hx.outs[0].phase == 'l'
    hx.V = None
    hx.T = 320
    hx.simulate()
    assert_allclose(hx.outs[0].T, 320, rtol=1e-3)
    hx.H = feed.H
    with pytest.raises(RuntimeError):
        hx.simulate() # Too many specifications
    hx.T = None
    hx.simulate()
    assert_allclose(hx.outs[0].H, feed.H, rtol=1e-2, atol=1)
    hx.H = None
    with pytest.raises(RuntimeError):
        hx.simulate() # No specifications

def test_hxutility_nonrigorous_specifications():
    bst.F.set_flowsheet('hx_nonrigorous')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, Ethanol=10, T=300)
    hx = bst.HXutility(ins=feed, T=320, H=1000)
    with pytest.raises(RuntimeError):
        hx.simulate()
    hx.T = None
    hx.simulate()
    assert_allclose(hx.outs[0].H, 1000)
    hx.H = None
    hx.T = 400
    hx.V = 1
    hx.simulate()
    assert hx.outs[0].phase == 'g'
    hx.V = 0.5
    with pytest.raises(RuntimeError):
        hx.simulate()
    hx.V = 1
    hx.T = 290
    with pytest.raises(ValueError):
        hx.simulate() # Cannot cool while boiling
    vapor = bst.Stream(Water=100, T=400, phase='g')
    hx = bst.HXutility(ins=vapor, T=410, V=0)
    with pytest.raises(ValueError):
        hx.simulate() # Cannot heat while condensing

def test_hxutility_heat_and_cool_only():
    bst.F.set_flowsheet('hx_only')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, T=350)
    hx = bst.HXutility(ins=feed, T=300, heat_only=True)
    hx.simulate()
    assert_allclose(hx.outs[0].T, feed.T)
    hx = bst.HXutility(ins=feed, T=400, cool_only=True)
    hx.simulate()
    assert_allclose(hx.outs[0].T, feed.T)
    empty = bst.Stream()
    hx = bst.HXutility(ins=empty, T=400)
    hx.simulate()
    assert hx.outs[0].isempty()

def test_hxutility_pressure_drop():
    bst.F.set_flowsheet('hx_dP')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, Ethanol=10, T=300, P=3 * 101325)
    hx = bst.HXutility(ins=feed, V=1, rigorous=True, estimate_pressure_drop=True)
    hx.simulate()
    assert hx.outs[0].P < feed.P
    assert hx.design_results['Pressure drop'] > 0
    # NOTE: estimate_pressure_drop with T specification raises UnboundLocalError
    # (vapor_fraction not defined in HXutility._run); not tested.
    hx = bst.HXutility(ins=feed, T=320, dP=1000)
    hx.simulate()
    assert_allclose(hx.outs[0].P, feed.P - 1000)

def test_hxutility_as_auxiliary_exchanger():
    bst.F.set_flowsheet('hx_auxiliary')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    a = bst.Stream(Water=100, T=300)
    b = bst.Stream(Ethanol=10, T=300)
    hx = bst.HXutility()
    hx.simulate_as_auxiliary_exchanger([a, b], duty=1e6, P_in=101325, P_out=101325)
    assert hx.outs[0].T > a.T
    assert_allclose(hx.net_duty, 1e6, rtol=0.1)
    hx = bst.HXutility()
    hx.simulate_as_auxiliary_exchanger([a, b], duty=1e6, vle=False, scale=2)
    assert hx.outs[0].T > a.T
    out = bst.Stream(Water=100, Ethanol=10, T=320)
    hx = bst.HXutility()
    hx.simulate_as_auxiliary_exchanger([a, b], [out], P_out=101325)
    assert hx.net_duty > 0
    hx = bst.HXutility()
    hx.simulate_as_auxiliary_exchanger([a, b], [out], duty=1e5)
    hx = bst.HXutility()
    hx.simulate_as_auxiliary_exchanger([a, b], [out], duty=1e5, vle=False)
    hx = bst.HXutility()
    with pytest.raises(ValueError):
        hx.simulate_as_auxiliary_exchanger([a, b])

def create_hxutilities(feed, T):
    # HXutilities does not define design inputs (e.g., T) like HXutility,
    # so specifications are set as attributes.
    hx = bst.HXutilities(ins=feed, rigorous=False)
    hx.T = T
    hx.V = hx.H = hx.dP = None
    return hx

def test_hxutilities():
    bst.F.set_flowsheet('hx_utilities')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, T=300)
    hx = create_hxutilities(feed, 450)
    hx.simulate()
    assert len(hx.heat_exchangers) >= 1
    assert_allclose(hx.outs[0].T, 450)
    feed = bst.Stream(Water=100, T=450, P=10*101325)
    hx = create_hxutilities(feed, 280)
    hx.simulate()
    assert len(hx.heat_exchangers) >= 1
    feed = bst.Stream(Water=100, T=300)
    hx = create_hxutilities(feed, 300)
    hx.simulate()
    assert hx.heat_exchangers == []

def test_hxprocess_misc():
    bst.F.set_flowsheet('hx_process')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    hot = bst.Stream(Water=100, T=360, P=2*101325)
    cold = bst.Stream(Ethanol=50, T=300, P=2*101325)
    hx = bst.HXprocess(ins=(hot, cold), dP0=1000, dP1=2000)
    hx.simulate()
    assert hx.design_results['Tube side pressure drop'] > 0
    assert hx.design_results['Shell side pressure drop'] > 0
    assert hx.Q > 0
    hx.Q = hx.Q
    assert len(hx.get_streams()) == 4
    assert_allclose(hx.outs[0].P, hot.P - 1000)
    # Empty inlets
    hx = bst.HXprocess(ins=(bst.Stream(), cold))
    hx.simulate()
    assert hx.total_heat_transfer == 0
    hx = bst.HXprocess(ins=(hot, bst.Stream()))
    hx.simulate()
    assert hx.total_heat_transfer == 0
    with pytest.raises(NotImplementedError):
        bst.HXprocess(ins=(hot, cold), estimate_pressure_drop=True)
    # Phenomena-oriented helpers
    hx = bst.HXprocess(ins=(hot, cold))
    hx.simulate()
    hx._update_nonlinearities()
    unit, C = hx._get_energy_departure_coefficient(hx.outs[0])
    assert unit is hx and C > 0
    T0 = hx.outs[0].T
    hx._update_energy_variable(1.)
    assert_allclose(hx.outs[0].T, T0 + 1.)

def test_hxprocess_reset_streams():
    bst.F.set_flowsheet('hx_process_reset')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    hot = bst.Stream(Water=100, T=360)
    cold = bst.Stream(Ethanol=50, T=300)
    M = bst.Mixer(ins=hot)
    hx = bst.HXprocess(ins=(M-0, cold))
    hx.reset_streams_at_setup = True
    sys = bst.System.from_units(units=[M, hx])
    sys.simulate()
    assert hx.total_heat_transfer > 0

def test_furnace():
    bst.F.set_flowsheet('furnace')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, T=300, phase='g')
    hx = bst.HXutility(ins=feed, T=700, rigorous=False)
    hx.simulate()
    assert hx.isfurnace
    assert 'Furnace' in hx.baseline_purchase_costs
    check_costs(hx)

if __name__ == '__main__':
    test_hxutility_exchanger_types()
    test_hxutility_rigorous_specifications()
    test_hxutility_nonrigorous_specifications()
    test_hxutility_heat_and_cool_only()
    test_hxutility_pressure_drop()
    test_hxutility_as_auxiliary_exchanger()
    test_hxutilities()
    test_hxprocess_misc()
    test_hxprocess_reset_streams()
    test_furnace()
