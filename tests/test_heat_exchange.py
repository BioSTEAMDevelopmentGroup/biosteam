# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests for heat exchanger units and the counter-current heat exchange solver.
"""
import biosteam as bst
from numpy.testing import assert_allclose
from biosteam.units.design_tools.heat_transfer import heat_exchange_to_condition

def test_heat_exchange_to_condition_respects_H_lim_at_bubble_point():
    # When the temperature limit coincides with the stream's bubble point
    # (within the solver's 1e-3 K tolerance), the outlet is set to the
    # saturated phase; the enthalpy limit must still be honored.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    s_in = bst.Stream('w_in', Water=2000., T=350., P=101325., phase='l',
                      units='kmol/hr')
    T_bp = s_in.bubble_point_at_P().T
    s_lim = s_in.copy(); s_lim.T = 360.
    H_lim = s_lim.H
    # heating: limit below full vaporization at the bubble point
    s_out = s_in.copy()
    Q = heat_exchange_to_condition(s_in, s_out, T=T_bp, H_lim=H_lim, heating=True)
    assert_allclose(s_out.H, H_lim, rtol=1e-9)
    assert_allclose(Q, H_lim - s_in.H, rtol=1e-9)
    assert s_out.phase == 'l' and s_out.T < T_bp
    # cooling: a saturated vapor whose limit is above full condensation
    v_in = s_in.copy('w_vap'); v_in.phase = 'g'; v_in.T = T_bp
    s_out = v_in.copy()
    H_lim = v_in.H - 0.5 * (v_in.H - s_in.H)   # half-way to liquid at 350 K
    Q = heat_exchange_to_condition(v_in, s_out, T=T_bp, H_lim=H_lim, heating=False)
    assert_allclose(s_out.H, H_lim, rtol=1e-9)
    assert_allclose(Q, H_lim - v_in.H, rtol=1e-9)

def test_HXprocess_H_lim_when_pinch_is_at_bubble_point():
    # Hot liquid water at exactly T_bp + dT against cold liquid water with an
    # enthalpy limit below vaporization: the cold outlet may not exceed its
    # enthalpy limit just because its temperature limit lands on the
    # bubble point.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    cold = bst.Stream('cold', Water=2000., T=350., P=101325., phase='l',
                      units='kmol/hr')
    T_bp = cold.bubble_point_at_P().T
    hot = bst.Stream('hot', Water=1000., T=T_bp + 5., P=5e5, phase='l',
                     units='kmol/hr')
    s_lim = cold.copy(); s_lim.T = 360.
    H_lim = s_lim.H
    hx = bst.HXprocess('hx', ins=(cold, hot), H_lim0=H_lim, T_lim1=355., dT=5.)
    hx.simulate()
    cold_out, hot_out = hx.outs
    assert_allclose(cold_out.H, H_lim, rtol=1e-9)
    assert_allclose(hx.Q, H_lim - cold.H, rtol=1e-9)
    assert_allclose(hot_out.H, hot.H - hx.Q, rtol=1e-9)
    assert hot_out.T > 355.   # the hot stream was not the limiting side

def test_HXprocess_never_modifies_inlets():
    # A superheated-liquid hot inlet (ethanol at 370 K, 1 atm; bp 351.4 K)
    # against a two-phase cold inlet of the same fluid: the solver finds no
    # feasible exchange. It must leave both inlet streams untouched rather
    # than overwrite one with its (already modified) outlet.
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    hot = bst.Stream('hot', Ethanol=700., T=370., P=101325., phase='l',
                     units='kmol/hr')
    cold = bst.Stream('cold', Ethanol=400., T=348.6, P=101325., phase='l',
                      units='kmol/hr')
    cold.vle(H=cold.H + 4e6, P=101325.)   # two-phase at the bubble point
    H_hot, H_cold = hot.H, cold.H
    hx = bst.HXprocess('hx', ins=(cold, hot), dT=5.)
    hx._run()
    assert hx.Q == 0.
    assert_allclose([hot.H, cold.H], [H_hot, H_cold], rtol=1e-12)
    assert hot.phase == 'l' and hot.T == 370.

if __name__ == '__main__':
    test_heat_exchange_to_condition_respects_H_lim_at_bubble_point()
    test_HXprocess_H_lim_when_pinch_is_at_bubble_point()
    test_HXprocess_never_modifies_inlets()
