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
from numpy.testing import assert_allclose
from thermo import SRK, PR


def test_isentropic_helium_turbine():
    bst.settings.set_thermo(["He"], cache=True)
    eta = 1.0
    feed = bst.Stream(He=1, T=-160 + 273.15, P=20e5, phase='g')
    K = bst.units.IsentropicTurbine(ins=feed, P=5e5, eta=eta)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, -208.16+273.15, 5e5],
        rtol=1e-4,
    )
    # check compressor design
    ideal_power = -0.278
    eta_motor = K.design_results['Driver efficiency']
    expected_power = ideal_power * eta * eta_motor
    actual_power = K.power_utility.rate
    assert_allclose(
        [actual_power, K.design_results['Ideal power'],
           K.design_results['Ideal duty'], K.design_results['Turbines in parallel']],
        [expected_power, ideal_power, 0, 1],
        rtol=1e-3,
    )
    pass


if __name__ == '__main__':
    test_isentropic_helium_turbine()
