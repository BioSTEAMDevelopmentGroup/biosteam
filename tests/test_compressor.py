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
from thermo import SRK, PR


def test_isentropic_hydrogen_compressor():
    bst.settings.set_thermo(["H2"])
    bst.settings.chemicals.H2.V.g.method_P = 'IDEAL'
    eta = 0.7
    feed = bst.Stream(H2=1, T=25 + 273.15, P=101325, phase='g')
    K = bst.units.IsentropicCompressor(ins=feed, P=50e5, eta=eta)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, 1151.9888331779612, 5000000.0],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Centrifugal"
    assert K.design_results["Driver"] == "Steam turbine"
    ideal_power = 4.920933056717339
    eta_motor = K.baseline_cost_algorithms[K.design_results["Type"]].efficiencies["Steam turbine"]
    expected_power = ideal_power / eta / eta_motor
    actual_power = K.power_utility.consumption
    assert_allclose(
        [actual_power, K.design_results['Ideal power'],
         K.design_results['Ideal duty'], K.design_results['Compressors in parallel']],
        [expected_power, ideal_power, 0, 1],
        rtol=1e-3,
    )
    pass


def test_isentropic_two_phase_steam_compressor():
    bst.settings.set_thermo(["H2O"])
    bst.settings.chemicals.H2O.V.g.method_P = 'IDEAL'
    # feed is steam-water-mixture
    feed = bst.MultiStream(T=372.75, P=1e5, l=[('H2O', 0.1)], g=[('H2O', 0.9)])
    # from compress until superheated steam (test case for vle parameter)
    eta = 1
    K = bst.units.IsentropicCompressor(ins=feed, P=100e5, eta=eta, vle=True)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, 798.6308938602593, 10000000.0],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Centrifugal"
    assert K.design_results["Driver"] == "Steam turbine"
    ideal_power = 5.4152214427108305
    eta_motor = K.baseline_cost_algorithms[K.design_results["Type"]].efficiencies["Steam turbine"]
    expected_power = ideal_power / eta / eta_motor
    actual_power = K.power_utility.consumption
    assert_allclose(
        [actual_power, K.design_results['Ideal power'],
           K.design_results['Ideal duty'], K.design_results['Compressors in parallel']],
        [expected_power, ideal_power, 0, 1],
        rtol=1e-3,
    )
    pass


def test_isothermal_hydrogen_compressor():
    thermo = bst.Thermo([bst.Chemical('H2', eos=SRK)])
    thermo.chemicals.H2.V.g.method_P = 'IDEAL'
    thermo.mixture.include_excess_energies = True
    bst.settings.set_thermo(thermo)
    feed = bst.Stream(H2=1, T=298.15, P=20e5, phase='g')

    ideal_power = 2.102711178815804
    ideal_duty = -7260.634742550824

    # eta = 1
    eta = 1
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5, eta=eta)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, feed.T, 350e5],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Reciprocating"
    assert K.design_results["Driver"] == "Electric motor"
    eta_motor = K.baseline_cost_algorithms[K.design_results["Type"]].efficiencies["Electric motor"]
    expected_power = ideal_power / eta / eta_motor
    actual_power = K.power_utility.rate
    assert_allclose(
        [actual_power, K.design_results['Ideal power'],
           K.design_results['Ideal duty'], K.design_results['Compressors in parallel']],
        [expected_power, ideal_power, ideal_duty, 1],
        rtol=1e-3,
    )
    # check heat utility
    expected_duty = ideal_duty / eta
    assert_allclose(
        [K.heat_utilities[0].unit_duty, K.heat_utilities[0].duty, K.heat_utilities[0].flow],
        [expected_duty, expected_duty, 7.5262772846646],
        rtol=1e-3,
    )

    # repeat with eta=0.7
    eta = 0.7
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5, eta=eta)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, feed.T, 350e5],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Reciprocating"
    assert K.design_results["Driver"] == "Electric motor"
    eta_motor = K.baseline_cost_algorithms[K.design_results["Type"]].efficiencies["Electric motor"]
    expected_power = ideal_power / eta / eta_motor
    actual_power = K.power_utility.rate
    assert_allclose(
        [actual_power, K.design_results['Ideal power'],
           K.design_results['Ideal duty'], K.design_results['Compressors in parallel']],
        [expected_power, ideal_power, ideal_duty, 1],
        rtol=1e-3,
    )
    # check heat utility
    expected_duty = ideal_duty / eta
    hu = K.heat_utilities[0]
    assert_allclose(
        [hu.unit_duty, hu.duty, hu.flow],
        [expected_duty, expected_duty, 7.526277284664599 / eta],
        rtol=1e-3,
    )
    pass

def test_polytropic_hydrogen_compressor():
    thermo = bst.Thermo([bst.Chemical('H2', eos=PR)])
    thermo.chemicals.H2.V.g.method_P = 'IDEAL'
    thermo.mixture.include_excess_energies = True
    bst.settings.set_thermo(thermo)
    feed = bst.Stream(H2=1, T=25 + 273.15, P=20e5, phase='g')
    P = 350e5
    # eta=1 outlet state should be identical to ideal isentropic compression
    # schultz method
    eta = 1.0
    K = bst.units.PolytropicCompressor(ins=feed, P=P, eta=eta, method="schultz")
    K.simulate()
    out_poly = K.outs[0]
    K_isen = bst.units.IsentropicCompressor(ins=feed, P=P, eta=eta)
    K_isen.simulate()
    out_isen = K_isen.outs[0]
    assert len(K.heat_utilities) == 0 and K.net_duty == 0.
    assert_allclose(
        [out_poly.T, out_poly.P, out_poly.H, K.power_utility.rate],
        [out_isen.T, out_isen.P, out_isen.H, K_isen.power_utility.rate],
        rtol=1e-3,
    )
    # eta=1, hundseid method
    K = bst.units.PolytropicCompressor(ins=feed, P=P, eta=eta, method="hundseid")
    K.simulate()
    out_poly = K.outs[0]
    assert len(K.heat_utilities) == 0 and K.net_duty == 0.
    assert_allclose(
        [out_poly.T, out_poly.P, out_poly.H, K.power_utility.rate],
        [out_isen.T, out_isen.P, out_isen.H, K_isen.power_utility.rate],
        rtol=1e-3,
    )
    # eta=0.7, hundseid method
    eta=0.7
    feed = bst.Stream(H2=1, T=25 + 273.15, P=20e5, phase='g')
    K = bst.units.PolytropicCompressor(ins=feed, P=P, eta=eta, method="hundseid", n_steps=200)
    K.simulate()
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, 958.0835733924658, P],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Reciprocating"
    assert K.design_results["Driver"] == "Electric motor"
    assert len(K.heat_utilities) == 0 and K.net_duty == 0.
    expected_power = 6.484199755124557
    actual_power = K.power_utility.rate
    assert_allclose(
        [actual_power, K.design_results['Compressors in parallel']],
        [expected_power, 1],
        rtol=1e-3,
    )
    pass


def test_multistage_hydrogen_compressor_simple():
    thermo = bst.Thermo([bst.Chemical('H2', eos=PR)])
    thermo.chemicals.H2.V.g.method_P = 'IDEAL'
    thermo.mixture.include_excess_energies = True
    bst.settings.set_thermo(thermo)
    feed = bst.Stream("feed", H2=1, T=25 + 273.15, P=20e5, phase='g')
    P = 350e5
    # test simple setup
    n_stages = 5
    pr = (P / feed.P) ** (1 / n_stages)
    K = bst.units.MultistageCompressor(ins=feed, outs=('outlet'), pr=pr, n_stages=n_stages, eta=0.7)
    K.simulate()
    # check inlet and outlet
    assert K.ins[0].ID == "feed"
    assert K.compressors[0].ins[0] == K.ins[0]
    assert K.outs[0].ID == "outlet"
    assert K.hxs[-1].outs[0] == K.outs[0]
    # check outlet state
    out = K.outs[0]
    assert_allclose(
        [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
        [1.0, 0.0, feed.T, P],
        rtol=1e-3,
    )
    # check compressor design
    assert K.design_results["Type"] == "Multistage compressor"
    # assert K.design_results["Driver"] == "Electric motor"
    assert_allclose(
        [K.design_results["Area"],
         K.design_results['Tube side pressure drop'],
         K.design_results['Shell side pressure drop']],
        [1.7906277579485619, 15.0, 25.0],
        rtol=1e-3,
    )
    # check heat utilities
    heat_utilities = bst.HeatUtility.sum_by_agent(K.heat_utilities)
    assert heat_utilities[1].ID == 'chilled_water'
    assert_allclose(
        [len(heat_utilities), heat_utilities[1].duty, heat_utilities[1].flow, heat_utilities[1].cost],
        [2, -11360.26358696248, 7.529257798470469, 0.056801317934812405],
        rtol=1e-3,
    )
    # check power utility
    assert_allclose(
        [K.power_utility.consumption, K.power_utility.production, K.power_utility.rate, K.power_utility.cost],
        [4.669790810710772, 3.86967926182525, 0.8001115488855222, 0.06256872312284784],
        rtol=1e-3,
    )
    pass


def test_multistage_hydrogen_compressor_advanced():
    thermo = bst.Thermo([bst.Chemical('H2', eos=PR)])
    thermo.chemicals.H2.V.g.method_P = 'IDEAL'
    thermo.mixture.include_excess_energies = True
    bst.settings.set_thermo(thermo)
    feed = bst.Stream("feed", H2=1, T=25 + 273.15, P=20e5, phase='g')
    P = 350e5
    # test advanced setup
    Ps = [60e5, 90e5, 200e5, 300e5, 350e5]
    Ts = [450, 400, 350, 300, feed.T]
    ks = [
        bst.units.IsentropicCompressor(P=Ps[0], eta=0.6),
        bst.units.PolytropicCompressor(P=Ps[1], eta=0.7),
        bst.units.IsentropicCompressor(P=Ps[2], eta=0.75),
        bst.units.PolytropicCompressor(P=Ps[3], eta=0.8),
        bst.units.IsentropicCompressor(P=Ps[4], eta=0.85),
    ]
    hxs = [bst.units.HXutility(T=T) for T in Ts]
    K = bst.units.MultistageCompressor(ins=feed, outs=('outlet'), compressors=ks, hxs=hxs)
    K.simulate()
    # check inlet and outlet
    assert K.ins[0].ID == "feed"
    assert K.compressors[0].ins[0] == K.ins[0]
    assert K.outs[0].ID == "outlet"
    assert K.hxs[-1].outs[0] == K.outs[0]
    # check hx outlet states
    for hx, P, T in zip(K.hxs, Ps, Ts):
        out = hx.outs[0]
        assert_allclose(
            [out.vapor_fraction, out.liquid_fraction, out.T, out.P],
            [1.0, 0.0, T, P],
            rtol=1e-3,
        )
    # check compressor design
    assert K.design_results["Type"] == "Multistage compressor"
    # assert K.design_results["Driver"] == "Electric motor"
    assert_allclose(
        [K.design_results["Area"],
         K.design_results['Tube side pressure drop'],
         K.design_results['Shell side pressure drop']],
        [1.1008174238336854, 15.0, 25.0],
        rtol=1e-3,
    )
    # check heat utilities
    heat_utilities = bst.HeatUtility.sum_by_agent(K.heat_utilities)
    assert len(heat_utilities) == 3
    assert heat_utilities[2].ID == 'chilled_water'
    assert heat_utilities[1].ID == 'cooling_water'
    assert heat_utilities[0].ID == 'high_pressure_steam'
    assert_allclose(
        [heat_utilities[2].duty, heat_utilities[2].flow, heat_utilities[2].cost],
        [-3694.721706001324, 2.4487558765814454, 0.018473608530006624],
        rtol=1e-3,
    )
    assert_allclose(
        [heat_utilities[1].duty, heat_utilities[1].flow, heat_utilities[1].cost],
        [-10377.599222102164, 7.087204175252751, 0.003457492556897055],
        rtol=1e-3,
    )
    # check power utility
    assert_allclose(
        [K.power_utility.consumption, K.power_utility.production, K.power_utility.rate, K.power_utility.cost],
        [6.021683407067501, 5.848476544603433, 0.1732068624640677, 0.013544776644690094],
        rtol=1e-3,
    )

def test_compressor_design():
    bst.settings.set_thermo(["H2"], cache=True)
    bst.settings.chemicals.H2.V.g.method_P = 'IDEAL'
    feed = bst.Stream(H2=1, T=298.15, P=1e5, phase='g')
    K = bst.units.IsothermalCompressor(ins=feed, P=350e5)
    
    # Test default compressor types across pressures:
    psig_to_pascal = lambda psig: (psig + 14.6959) * 101325 / 14.6959
    def compressor_type_at_psig(psig):
        K.P = psig_to_pascal(psig)
        K.simulate()
        return K.design_results["Type"]
    
    # Pressures [psig] in the same order of magnitude as upper bound for each
    # compressor type.
    pressures = [1e2, 1e3, 1e4] 
    compressor_types = ['Screw', 'Centrifugal', 'Reciprocating']
    assert [compressor_type_at_psig(i) for i in pressures] == compressor_types
    
    # Test compressors in parallel:
    # Actual cubic feet per minute in the same order of magnitude as upper
    # bound for each compressory type.
    acfms = [1.5e4, 1e5, 5e3]
    
    def get_compressors_in_parallel(psig, acfm):
        K.ins[0].set_total_flow(acfm, 'cfm')
        K.P = psig_to_pascal(psig)
        K.simulate()
        return K.design_results['Compressors in parallel']
    
    compressors_in_parallel = [1 , 2]
    for psig, acfm in zip(pressures, acfms):
        assert compressors_in_parallel == [get_compressors_in_parallel(psig, i * acfm) for i in compressors_in_parallel]
        
    # Test user overwrite:
    K.P = psig_to_pascal(1e3) 
    K.compressor_type = "Centrifugal"
    K.driver = 'Electric motor'
    K.driver_efficiency = 0.9
    K.simulate()
    assert K.design_results['Type'] == 'Centrifugal'
    assert K.design_results['Driver'] == 'Electric motor'
    assert K.design_results['Driver efficiency'] == 0.9
    
    # Test capital cost correlations and factors
    def installed_equipment_cost_at_psig(psig):
        K.P = psig_to_pascal(psig)
        K.simulate()
        return K.installed_cost
    
    K.driver = K.driver_efficiency = K.compressor_type = 'Default'
    assert_allclose([installed_equipment_cost_at_psig(i) for i in pressures],
                    [ 2267264.256152,  5179043.514644, 10479970.002024],
                    rtol=1e-3)

def test_multistage_setup_does_not_recreate_subcomponents():
    bst.settings.set_thermo(["H2"])
    bst.settings.chemicals.H2.V.g.method_P = 'IDEAL'
    feed = bst.Stream("feed", H2=1, T=25 + 273.15, P=20e5, phase='g')
    P = 350e5
    n_stages = 5
    pr = (P / feed.P) ** (1 / n_stages)
    K = bst.units.MultistageCompressor(ins=feed, outs=('outlet'), pr=pr, n_stages=n_stages, eta=0.7)
    K._setup()
    a = K.compressors[0]
    K._setup()
    b = K.compressors[0]
    assert a==b
    
def test_multistage_setup_updates_after_changing_specifications():
    bst.settings.set_thermo(["H2"])
    bst.settings.chemicals.H2.V.g.method_P = 'IDEAL'
    feed = bst.Stream("feed", H2=1, T=25 + 273.15, P=20e5, phase='g')
    P = 350e5
    n_stages = 5
    pr = (P / feed.P) ** (1 / n_stages)
    K = bst.units.MultistageCompressor(ins=feed, outs=('outlet'), pr=pr, n_stages=n_stages, eta=0.7)
    K._setup()
    units_5 = (K.compressors, K.hxs)
    for i in units_5: assert len(i) == 5
    K.n_stages = 6
    K._setup()
    units_6 = (K.compressors, K.hxs)
    for i in units_6: assert len(i) == 6


if __name__ == '__main__':
    test_compressor_design()
    test_isentropic_hydrogen_compressor()
    test_isentropic_two_phase_steam_compressor()
    test_isothermal_hydrogen_compressor()
    test_polytropic_hydrogen_compressor()
    test_multistage_hydrogen_compressor_simple()
    test_multistage_hydrogen_compressor_advanced()
    test_multistage_setup_does_not_recreate_subcomponents()
    test_multistage_setup_updates_after_changing_specifications()