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


def test_isenthalpic_methane_liquefaction():
    thermo = bst.Thermo([bst.Chemical('CH4')])
    thermo.mixture.include_excess_energies = True
    bst.settings.set_thermo(thermo, cache=True)
    feed = bst.Stream(CH4=1, T=-63 + 273.15, P=100e5, phase='g')
    V = bst.units.IsenthalpicValve(ins=feed, P=1e5, vle=True)
    V.simulate()
    # check outlet state
    out = V.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P, out.H],
        [0.8353, 0.1646, 111.474, 1e5, feed.H],
        rtol=1e-2,
    )
    pass


if __name__ == '__main__':
    test_isenthalpic_methane_liquefaction()
