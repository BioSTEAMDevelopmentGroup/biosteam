# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biosteam.units import design_tools as d
from numpy.testing import assert_allclose

def reset_CEPCI(f):
    def g():
        try:
            f()
        finally:
            bst.default_CEPCI()
    g.__name__ = f.__name__
    return g

@reset_CEPCI
def test_electric_motor_cost():
    assert_allclose(
        d.electric_motor_cost(10),
        620.5362044090548, 
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.electric_motor_cost(10),
        874.7646934400772, 
    )

@reset_CEPCI
def test_horizontal_vessel_purchase_cost():
    assert_allclose(
        d.compute_horizontal_vessel_purchase_cost(W=1e3),
        8857.95780922862, 
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_horizontal_vessel_purchase_cost(W=1e3),
        12486.988982172503,
    )

@reset_CEPCI
def test_vertical_vessel_platform_and_ladders_purchase_cost():
    assert_allclose(
        d.compute_vertical_vessel_platform_and_ladders_purchase_cost(3, 10),
        4708.532188390147, 
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_vertical_vessel_platform_and_ladders_purchase_cost(3, 10),
        6637.578415351749,
    )

@reset_CEPCI
def test_vertical_vessel_purchase_cost():
    assert_allclose(
        d.compute_vertical_vessel_purchase_cost(1e3),
        13319.089264026192,
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_vertical_vessel_purchase_cost(1e3),
        18775.808654133838,
    )

@reset_CEPCI
def test_purchase_cost_of_trays():
    assert_allclose(
        d.compute_purchase_cost_of_trays(10, 6),
        17092.681096675424,
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_purchase_cost_of_trays(10, 6),
        24095.409475489585,
    )

@reset_CEPCI
def test_empty_tower_cost():
    assert_allclose(
        d.compute_empty_tower_cost(3000),
        30829.15740696664,
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_empty_tower_cost(3000),
        43459.60515519526,
    )

@reset_CEPCI
def test_plaform_ladder_cost():
    assert_allclose(
        d.compute_plaform_ladder_cost(6, 30),
        16225.103055811393,
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_plaform_ladder_cost(6, 30),
        22872.39197294998,
    )

@reset_CEPCI
def test_closed_vessel_turbine_purchase_cost():
    assert_allclose(
        d.compute_closed_vessel_turbine_purchase_cost(power=10),
        15264.970467626266,
    )
    bst.settings.CEPCI = 800
    assert_allclose(
        d.compute_closed_vessel_turbine_purchase_cost(power=10),
        21518.901099737468,
    )
    
if __name__ == '__main__':
    test_electric_motor_cost()
    test_horizontal_vessel_purchase_cost()
    test_vertical_vessel_platform_and_ladders_purchase_cost
    test_vertical_vessel_purchase_cost()
    test_purchase_cost_of_trays()
    test_empty_tower_cost()
    test_plaform_ladder_cost()
    test_closed_vessel_turbine_purchase_cost()