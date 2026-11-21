# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import matplotlib
matplotlib.use('Agg')
import pytest
import biosteam as bst
import numpy as np
from numpy.testing import assert_allclose
from biosteam.units.design_tools import (
    aeration,
    flash_vessel_design as fvd,
    heat_transfer as ht,
    vacuum,
    Gibbs_equilibrium_reaction as ger,
)
from biosteam.units.design_tools.pressure_vessel import PressureVessel

def test_flash_vessel_design_tables():
    for DRho in (5, 30, 60):
        for Hlr in (12, 24, 30):
            G = fvd.GTable(DRho, Hlr)
            assert G > 0 and np.isfinite(G)
    assert_allclose(fvd.HNATable(1, 0.5), 0.5, rtol=1e-2)
    assert_allclose(fvd.HNATable(2, 0.5), 0.5, rtol=1e-2)
    assert fvd.compute_low_liq_level_height(1, 100, 5) == 1.25
    assert fvd.compute_low_liq_level_height(1, 400, 5) == 0.5
    heights = [fvd.compute_low_liq_level_height(2, 100, D) for D in (3, 5, 8, 10, 12, 20)]
    assert heights == sorted(heights)
    assert_allclose(heights[0], 0.75)
    assert_allclose(heights[-1], 1.25)

def test_aeration_diffusion_coefficients():
    assert aeration.LeBas_volume({'C': 1, 'O': 2}) == 14.8 + 2 * 7.4
    with pytest.raises(KeyError):
        aeration.LeBas_volume({'Xx': 1})
    assert aeration.Vm_at_Tb('O2') == 25.6
    V_glucose = aeration.Vm_at_Tb('Glucose')
    assert_allclose(V_glucose, aeration.LeBas_volume({'C': 6, 'H': 12, 'O': 6}))
    D = aeration.HaydukLaudie_diffusion_coefficient(T=298.15)
    assert 1e-5 < D < 1e-4 # cm2/s
    D_CO2 = aeration.HaydukLaudie_diffusion_coefficient(T=298.15, P=101325, A='CO2')
    assert D_CO2 < D # CO2 is larger than O2
    D = aeration.WilkeChang_diffusion_coefficient(T=298.15, A='O2')
    assert 1e-5 < D < 1e-4 # cm2/s
    D = aeration.WilkeChang_diffusion_coefficient(T=298.15, A='Ethanol', B='Water')
    assert D > 0 and np.isfinite(D)
    D = aeration.WilkeChang_diffusion_coefficient(T=298.15, V_A=25.6, mu_B=0.89, phi_B=2.6, M_B=18.)
    assert 1e-5 < D < 1e-4 # cm2/s

def test_aeration_kLa_correlations():
    with pytest.raises(ValueError):
        aeration.register(lambda: None) # Invalid name
    kLa = aeration.kLa_stirred_Riet(P=1e4, V=10, U=0.05)
    kLa_Riet = aeration.kLa_stirred_Riet(P=1e4, V=10, U=0.05, coefficients="Van't Riet")
    kLa_custom = aeration.kLa_stirred_Riet(P=1e4, V=10, U=0.05, coefficients=(0.026, 0.4, 0.5))
    assert kLa > 0 and kLa_Riet > 0
    assert_allclose(kLa_Riet, kLa_custom)
    assert_allclose(aeration.P_at_kLa_Riet(kLa, V=10, U=0.05), 1e4, rtol=1e-6)
    assert_allclose(aeration.P_at_kLa_Riet(kLa_Riet, V=10, U=0.05, coefficients="Van't Riet"), 1e4, rtol=1e-6)
    assert_allclose(aeration.P_at_kLa_Riet(kLa_Riet, V=10, U=0.05, coefficients=(0.026, 0.4, 0.5)), 1e4, rtol=1e-6)
    assert aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1, system_like='Glycerol') > 0
    assert aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1, system_like='Xanthan') > 0
    assert aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1, system_like='Newtonian', coefficients=(1e-2, 0.5, -0.1)) > 0
    assert aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1, system_like='Non-Newtonian', coefficients=(5e-3, 0.6, -0.5, -0.8)) > 0
    with pytest.raises(ValueError):
        aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1)
    with pytest.raises(ValueError):
        aeration.kla_bubcol_DeJesus(Q=1, mu=0.01, k=1, n=1, system_like='Other', coefficients=(1, 1, 1))
    kwargs = dict(D=0.15, mu_l=1e-3, rho_l=1000, D_l=2e-9, g=9.81, sigma_l=0.07, epsilon_g=0.1)
    assert aeration.kla_bubcol_Akita_Yoshida(**kwargs) > 0
    assert aeration.kla_bubcol_Posarac_Tekic(**kwargs) > 0
    assert aeration.kla_bubcol_Seno(**kwargs, V_g=0.02, V_l=0.01) > 0
    assert aeration.kla_bubcol_Suh(**kwargs, V_g=0.02, V_l=0.01) > 0
    assert aeration.kla_bubcol_Dewes(V_g=0.02, mu_l=1e-3, rho_g=1.2) > 0
    assert aeration.kla_bubcol_Shah(D=0.29, V_l=0.001, g=9.81, V_g=0.05, rho_l=1000, mu_l=1e-3) > 0
    assert aeration.kla_stirred_Labik(N=5, D=0.15, v_s=5e-3, P0=1000) > 0
    organisms = ('P.shermanii', 'S.cerevisiae', 'P.chrysogenum_pellet', 'P.chrysogenum_mycelia')
    for organism in organisms:
        kLa = aeration.kLa_stirred_Galaction(
            'Submerged', V=1, C_x=10, v_s=0.01, P_a=1000, organism_type=organism
        )
        assert kLa > 0 and np.isfinite(kLa)
        kLa = aeration.kLa_stirred_Galaction(
            'Surface', V=1, C_x=10, v_s=0.01, P=1000, organism_type=organism
        )
        assert kLa > 0 and np.isfinite(kLa)
    for aeration_type in ('Submerged', 'Surface'):
        with pytest.raises(ValueError):
            aeration.kLa_stirred_Galaction(aeration_type, V=1, C_x=10, v_s=0.01, P_a=1, P=1, organism_type='E.coli')
    with pytest.raises(ValueError):
        aeration.kLa_stirred_Galaction('Other', V=1, C_x=10, v_s=0.01)
    assert aeration.log_mean_driving_force(0.2, 0.3, 0.2) > 0
    assert aeration.log_mean_driving_force(0.3, 0.4, 0.1, 0.05) > 0

def test_vacuum_design_tools():
    systems = vacuum.get_preferred_vacuum_systems('Dry-vacuum pump')
    assert list(systems)[0] == 'Dry-vacuum pump'
    systems = vacuum.get_preferred_vacuum_systems(['Steam-jet ejector'])
    assert list(systems) == ['Steam-jet ejector']
    with pytest.raises(ValueError):
        vacuum.get_preferred_vacuum_systems('Not a vacuum')
    with pytest.raises(ValueError):
        vacuum.get_preferred_vacuum_systems(['Not a vacuum'])
    available = vacuum.get_available_vacuum_systems(100, 100)
    assert len(available) > 0
    defaults = vacuum._default_vacuum_systems
    name, grade, N = vacuum.select_vacuum_system(defaults, 1e7, 100)
    assert N > 1
    name, grade, N = vacuum.select_vacuum_system(defaults, 1e-3, 100, ignore_F_lb=True)
    assert N == 1
    with pytest.raises(bst.exceptions.DesignError):
        vacuum.select_vacuum_system(defaults, 100, 0)
    for name, vacuum_sys in defaults.items():
        if name == 'Rotary-vane pump': continue
        for grade in vacuum_sys:
            cost = vacuum.calculate_vacuum_cost(name, grade, 100, 100, 100)
            assert cost > 0
    results = vacuum.compute_vacuum_system_power_and_cost(
        F_mass=100, F_vol=100, P_suction=5000, vessel_volume=0,
    )
    assert results
    for preference in ('Liquid-ring pump', 'Dry-vacuum pump', 'Rotary-vane pump'):
        results = vacuum.compute_vacuum_system_power_and_cost(
            F_mass=1, F_vol=1, P_suction=5000, vessel_volume=1,
            vacuum_system_preference=preference,
        )
        assert results

def test_heat_transfer_tools():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    assert ht.heuristic_pressure_drop(0, 1) == 1.5
    assert ht.heuristic_pressure_drop(0, 0) == 5.
    assert ht.heuristic_pressure_drop(1, 1) == 3.
    ci = bst.Stream(Water=1, T=300)
    co = bst.Stream(Water=1, T=320)
    hi = bst.Stream(Water=1, T=400, phase='g')
    ho = bst.Stream(Water=1, T=370)
    assert ht.heuristic_tubeside_and_shellside_pressure_drops(ci, hi, co, ho) == (1.5, 5.)
    assert ht.heuristic_tubeside_and_shellside_pressure_drops(ci, hi, co, ho, False) == (5., 1.5)
    # Heat exchange to condition with enthalpy limits
    s_in = bst.Stream(Water=1, T=300)
    s_out = s_in.copy()
    assert ht.heat_exchange_to_condition(s_in, s_out, T=350, H_lim=s_in.H - 1, heating=True) == 0.
    assert ht.heat_exchange_to_condition(s_in, s_out, T=280, H_lim=s_in.H + 1, heating=False) == 0.
    s_out = s_in.copy()
    H_lim = s_in.H + 100
    Q = ht.heat_exchange_to_condition(s_in, s_out, T=350, phase='l', H_lim=H_lim, heating=True)
    assert_allclose(Q, 100, rtol=1e-3)
    s_in = bst.Stream(Water=1, T=350)
    s_out = s_in.copy()
    H_lim = s_in.H - 100
    Q = ht.heat_exchange_to_condition(s_in, s_out, T=300, phase='l', H_lim=H_lim, heating=False)
    assert_allclose(Q, -100, rtol=1e-3)

def test_counter_current_heat_exchange_limits():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    hot = bst.Stream(Water=1, Ethanol=1, T=360)
    cold = bst.Stream(Water=1, Ethanol=1, T=300)
    hot_out = hot.copy()
    cold_out = cold.copy()
    # Temperature limits beyond feasibility give no heat exchange
    Q = ht.counter_current_heat_exchange(hot, cold, hot_out, cold_out, 5, T_lim0=370)
    assert Q == 0.
    Q = ht.counter_current_heat_exchange(hot, cold, hot_out, cold_out, 5, T_lim1=290)
    assert Q == 0.
    Q = ht.counter_current_heat_exchange(hot, cold, hot_out, cold_out, 5, T_lim0=340, T_lim1=320)
    assert Q > 0
    assert hot_out.T >= 340 - 1e-3
    assert cold_out.T <= 320 + 1e-3
    # Too close for heat exchange
    Q = ht.counter_current_heat_exchange(hot, cold, hot_out, cold_out, 100)
    assert Q == 0.

class SimplePressureVessel(bst.Unit, PressureVessel):
    _N_outs = 1
    _units = PressureVessel._units

    def _default_vessel_type(self):
        return 'Horizontal'

    def _design(self):
        self._design_and_cost_pressure_vessel(
            1e6, 2, 6, pressure_units='Pa', length_units='m', n_vessels=2,
        )

def test_pressure_vessel():
    bst.F.set_flowsheet('pressure_vessel')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    V = SimplePressureVessel(ins=bst.Stream(Water=1))
    V.simulate()
    assert V.vessel_type == 'Horizontal'
    assert all([i > 0 for i in V.baseline_purchase_costs.values()])
    assert any(['x2' in i for i in V.baseline_purchase_costs])
    with pytest.warns(bst.exceptions.DesignWarning):
        V._horizontal_vessel_design(10, 5, 20)
    with pytest.warns(bst.exceptions.DesignWarning):
        V._vertical_vessel_design(10, 5, 20)
    design = V._vertical_vessel_design(100, 5, 20, annular_diameter=0.5, check=False)
    assert design['Jacketed diameter'] > 5.5
    design = V._horizontal_vessel_design(100, 5, 20, annular_diameter=0.5, check=False)
    assert design['Jacketed diameter'] > 5.5
    V._vessel_type = 'Spherical'
    with pytest.raises(RuntimeError):
        V._vessel_design(100, 5, 20)
    with pytest.raises(RuntimeError):
        V._vessel_purchase_cost(1e4, 5, 20)
    with pytest.raises(NotImplementedError):
        PressureVessel._default_vessel_type(V)

def test_minimize_Gibbs_free_energy():
    bst.settings.set_thermo(['N2O4', 'NO2'], pkg='ideal gas')
    feed = bst.Stream(N2O4=100, T=300, P=101325, phase='g')
    product = feed.copy()
    ger.minimize_Gibbs_free_energy(product, method='COBYLA')
    assert_allclose(product.F_mass, feed.F_mass, rtol=1e-3)
    assert product.imol['NO2'] > 0
    product = feed.copy()
    ger.minimize_Gibbs_free_energy(product, IDs=['N2O4', 'NO2'], method='COBYLA')
    assert_allclose(product.F_mass, feed.F_mass, rtol=1e-3)
    with pytest.raises(ValueError):
        ger.minimize_Gibbs_free_energy(feed.copy(), method='SHGO')
    import matplotlib.pyplot as plt
    ger.plot_Gibbs_equilibrium_reaction_surface()
    plt.close('all')

if __name__ == '__main__':
    test_flash_vessel_design_tables()
    test_aeration_diffusion_coefficients()
    test_aeration_kLa_correlations()
    test_vacuum_design_tools()
    test_heat_transfer_tools()
    test_counter_current_heat_exchange_limits()
    test_pressure_vessel()
    test_minimize_Gibbs_free_energy()
