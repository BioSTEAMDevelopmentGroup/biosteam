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

hr = bst.wastewater.high_rate

def set_wwt_thermo(flowsheet):
    bst.F.set_flowsheet(flowsheet)
    chemicals = bst.Chemicals(['CO2', 'AceticAcid', 'H2O', 'O2', 'N2', 'CH4'])
    bst.settings.set_thermo(hr.append_wwt_chemicals(chemicals), cache=True)

def test_internal_circulation_reactor():
    set_wwt_thermo('wwt_IC')
    inf = bst.Stream('inf', AceticAcid=1, H2O=100, units='kg/hr')
    R = hr.InternalCirculationRx('R1', ins=inf, outs=('biogas', 'eff', 'sludge'))
    R.simulate()
    # Reactions are not atomically balanced; check mass balance loosely
    assert abs(R.mass_balance_error()) < 0.01 * inf.F_mass
    assert R.installed_cost > 0
    assert 0 < R.organic_rm <= 1
    assert R.Se >= 0 and R.Sw >= 0 and R.Xw >= 0 and R.tau > 0
    assert R.i_rm is not None and R.decay_rxn is not None
    assert R.biogas_rxns and R.growth_rxns
    # "separate"-only parameters are None in "lumped" method
    for name in ('mu_max', 'b', 'Fxb', 'Fxt', 'Vb', 'Vt'):
        assert getattr(R, name) is None
    # Setters with validation
    R.q_Qw = 0.01
    R.q_Xw = 1.5
    R.mu_max = 0.01
    R.b = 0.001
    R.Fxb = 0.5
    R.Fxt = 0.5
    R.biodegradability = 0.8
    for name, value in [('method', 'other'), ('OLRall', -1), ('q_Qw', 2),
                        ('q_Xw', 0.5), ('mu_max', -1), ('b', -1),
                        ('Fxb', 2), ('Fxt', 2), ('biodegradability', 2)]:
        with pytest.raises(ValueError):
            setattr(R, name, value)
    with pytest.raises(ValueError):
        R.biodegradability = {'AceticAcid': 2}

def test_anaerobic_membrane_bioreactor():
    set_wwt_thermo('wwt_AnMBR')
    inf = bst.Stream('inf', AceticAcid=10, H2O=1000, units='kg/hr')
    R = hr.AnMBR('R2', ins=(inf, '', 'naocl', 'citric', 'bisulfite', 'air'),
                 outs=('biogas', 'perm', 'sludge', 'vent'),
                 include_pump_building_cost=True, include_excavation_cost=True)
    R.simulate()
    assert abs(R.mass_balance_error()) < 0.01 * inf.F_mass
    assert R.installed_cost > 0 and np.isfinite(R.installed_cost)
    assert 0 < R.organic_rm <= 1
    assert R.AF is None and R.AeF is None
    assert R.L_well == R.W_well == R.D_well == 0
    assert R.TMP_anaerobic > 0 # Note: TMP_aerobic getter is broken
    R.TMP_aerobic = 1.
    assert R.v_cross_flow > 0 and R.v_GAC == 0
    assert R.Q_gpm > 0 and R.Q_cfs > 0 and R.OLR > 0
    assert R.N_blower == 0 and R.N_degasser >= 1
    assert R.biogas_rxns is not None and R.growth_rxns is not None
    assert R.i_rm is not None and R.split is not None
    # Setters (assign the current values back)
    # Note: mod_per_cas_range and cas_per_tank_range setters are broken
    for name in ('N_train', 'cas_per_tank_spare', 'mod_surface_area', 'W_tank', 'D_tank',
                 'W_dist', 'W_eff', 'excav_slope', 'constr_access', 'HRT',
                 'recir_ratio', 'J_max', 'TMP_anaerobic', 'SGD',
                 'AFF', 'v_cross_flow', 'v_GAC', 'sludge_conc'):
        value = getattr(R, name)
        setattr(R, name, value)
    R.L_well = R.W_well = R.D_well = 8
    R.biodegradability = 0.9
    R.membrane_configuration = 'crossflow'
    R.membrane_type = 'multitube'
    R.membrane_material = 'pes'
    for name, value in [('reactor_type', 'other'), ('membrane_configuration', 'other'),
                        ('membrane_type', 'other'), ('membrane_material', 'other'),
                        ('biodegradability', 2)]:
        with pytest.raises(ValueError):
            setattr(R, name, value)
    with pytest.raises(ValueError):
        R.biodegradability = {'AceticAcid': 2}
    # Invalid designs
    for kwargs in [dict(include_aerobic_filter=True),
                   dict(membrane_type='hollow fiber'),
                   dict(add_GAC=True),
                   dict(membrane_material='sintered steel'),
                   dict(membrane_configuration='submerged', membrane_type='multi-tube'),
                   dict(membrane_configuration='submerged', membrane_type='hollow fiber',
                        membrane_material='ceramic'),
                   dict(membrane_configuration='submerged', membrane_type='flat sheet',
                        membrane_material='ceramic')]:
        with pytest.raises(bst.exceptions.DesignError):
            hr.AnMBR(None, ins=(bst.Stream(None, AceticAcid=10, H2O=1000, units='kg/hr'), '', '', '', '', ''), 
                     outs=('', '', '', ''), **kwargs)

def test_polishing_filter_properties():
    set_wwt_thermo('wwt_polishing_filter')
    inf = bst.Stream('inf', AceticAcid=1, H2O=100, units='kg/hr')
    PF = hr.PolishingFilter('PF', ins=(inf, '', 'air'), outs=('biogas', 'eff', 'waste', 'air_out'),
                            filter_type='aerobic')
    PF.simulate()
    assert abs(PF.mass_balance_error()) < 1e-6
    assert PF.installed_cost > 0
    assert 0 < PF.organic_rm <= 1
    assert PF.Q_cmd > 0 and PF.i_rm is not None and PF.split is not None
    for name in ('d_max', 'excav_slope', 'constr_access', 't_wall', 't_slab',
                 'recir_ratio', 'sludge_conc'):
        setattr(PF, name, getattr(PF, name))
    for name, value in [('OLR', -1), ('X_decomp', 2), ('X_growth', 2)]:
        with pytest.raises(ValueError):
            setattr(PF, name, value)

def test_sludge_handling():
    set_wwt_thermo('wwt_sludge_handling')
    feed = bst.Stream('sludge', WWTsludge=1, H2O=100, units='kg/hr')
    BT = hr.BeltThickener('BT', ins=feed, outs=('eff', 'sludge_out'), solubles=('H2O', 'AceticAcid'))
    BT.simulate()
    assert abs(BT.mass_balance_error()) < 1e-6
    # Too dry influent skips simulation
    feed.imass['H2O'] = 0.01
    with pytest.warns(UserWarning):
        BT.simulate()
    assert BT.outs[0].isempty()
    # Empty influent
    feed.empty()
    BT.simulate()
    assert BT.outs[1].isempty()

def test_system_units():
    set_wwt_thermo('wwt_system_units')
    from biosteam.wastewater.high_rate.system import (
        ReverseOsmosis, Skipped, CHP, BiogasUpgrading
    )
    RO = ReverseOsmosis('RO', ins=bst.Stream('ro_in', H2O=1000, AceticAcid=1, units='kg/hr'))
    RO.simulate()
    assert abs(RO.mass_balance_error()) < 1e-6
    assert RO.RO_treated_water.F_mass > RO.brine.F_mass > 0
    unit = bst.Pump('P_skipped', ins=bst.Stream('p_in', H2O=10, price=1.))
    unit.simulate()
    S = Skipped('S', ins=bst.Stream('s_in', H2O=10), outs=('s_out',),
                wwt_units=[unit], wwt_streams=[bst.Stream('s_extra', H2O=1)], clear_wwt=True)
    S.simulate()
    assert S.utility_cost == 0 and not unit.baseline_purchase_costs
    biogas = bst.Stream('biogas_in', CH4=10, CO2=5, units='kg/hr', phase='g')
    # Note: CHP is not tested; `get_combustion_energy` fails when combustion 
    # chemicals (e.g., P4O10, Ash) are missing.
    U = BiogasUpgrading('U', ins=(biogas, 'foo'), outs=('RNG', 'remained'), ratio=0.5)
    U.simulate()
    assert U.outs[0].imass['CH4'] > 0
    assert_allclose(U.outs[0].imass['CH4'], 10 * 0.5 * 0.95, rtol=1e-2)

def test_wwt_utils():
    set_wwt_thermo('wwt_utils')
    from biosteam.wastewater.high_rate import utils
    s = bst.Stream('s', AceticAcid=1, H2O=100, WWTsludge=1, units='kg/hr')
    assert utils.compute_stream_COD(bst.Stream()) == 0
    assert utils.compute_stream_COD(s) > 0
    assert utils.kph_to_tpd(s) > 0
    solubles = utils.get_soluble_IDs(s.chemicals, utils.default_insolubles)
    assert 'AceticAcid' in solubles
    with pytest.raises(ValueError):
        utils.get_digestion_rxns(s, 1., 0.9, 0.2, 'WWTsludge')
    assert utils.get_digestion_rxns(s, 0.5, 0.5, 0.1, 'WWTsludge')
    assert utils.get_BD_dct(s.chemicals, Extra=0.5)['Extra'] == 0.5
    assert 'WWTsludge' in utils.get_insoluble_IDs(s.chemicals, ['WWTsludge', 'Other'])
    class Dummy:
        Q_mgd = 100.
        recir_ratio = 1.
    pumps, building = utils.cost_pump(Dummy())
    assert pumps > 0 and building > 0

if __name__ == '__main__':
    test_internal_circulation_reactor()
    test_anaerobic_membrane_bioreactor()
    test_polishing_filter_properties()
    test_sludge_handling()
    test_system_units()
    test_wwt_utils()
