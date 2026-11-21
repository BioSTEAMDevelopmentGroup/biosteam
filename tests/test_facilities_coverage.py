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

def create_combustion_chemicals():
    chemicals = bst.Chemicals([
        'Water', 'CH4', 'O2', 'N2', 'CO2', 'SO2', 'H2S', 'CaO', 'CaSO4',
        bst.Chemical('Ash', search_db=False, default=True, phase='s', MW=1),
    ])
    for i in ('CaO', 'CaSO4'): chemicals[i].at_state('s')
    return chemicals

def create_boiler_system(cls, flowsheet, **kwargs):
    bst.settings.set_thermo(create_combustion_chemicals(), cache=True)
    bst.F.set_flowsheet(flowsheet)
    feed = bst.Stream('feed', Water=100, T=300)
    side_steam = bst.Stream('side_steam', Water=1, phase='g', T=400)
    with bst.System('sys') as sys:
        H = bst.HXutility('H', feed, T=360)
        B = cls('B', ins=[bst.Stream('solids', H2S=0.1), bst.Stream('gas')], 
                side_steam=side_steam, **kwargs)
    return sys, H, B

def check_boiler_results(sys, H, B):
    sys.simulate()
    # Steam from boiler satisfies the heater
    steam_duty = sum([i.duty for i in B.heat_utilities])
    assert_allclose(-steam_duty, H.net_duty, rtol=1e-2)
    assert B.installed_cost > 0 and np.isfinite(B.installed_cost)
    assert B.emissions.F_mass > 0
    assert B.ash_disposal.F_mass > 0
    assert B.makeup_water.F_mass > 0 and B.blowdown_water.F_mass > 0
    assert B.air.F_mass > 0
    assert not B.ins[4].isempty() # Lime

def test_boiler():
    sys, H, B = create_boiler_system(
        bst.Boiler, 'boiler_coverage',
        boiler_efficiency_basis='HHV', fuel_price=0.2, ash_disposal_price=0.02,
    )
    assert B.fuel_price == 0.2 and B.ash_disposal_price == 0.02
    check_boiler_results(sys, H, B)
    B.boiler_efficiency_basis = 'invalid'
    with pytest.raises(ValueError):
        sys.simulate()
    with pytest.raises(ValueError):
        bst.Boiler('B2', T_emissions=400)
    bst.Boiler('B3', natural_gas_price=0.2)

def test_boiler_turbogenerator_simple():
    sys, H, B = create_boiler_system(
        bst.BoilerTurbogenerator, 'boiler_turbogenerator_coverage',
        boiler_efficiency_basis='HHV', fuel_price=0.2, ash_disposal_price=0.02,
    )
    assert B.fuel_price == 0.2 and B.ash_disposal_price == 0.02
    check_boiler_results(sys, H, B)
    B.boiler_efficiency_basis = 'invalid'
    with pytest.raises(ValueError):
        sys.simulate()
    with pytest.raises(ValueError):
        bst.BoilerTurbogenerator('BT2', T_emissions=400)
    bst.BoilerTurbogenerator('BT3', natural_gas_price=0.2)

def test_boiler_without_lime():
    bst.settings.set_thermo(['Water', 'CH4', 'O2', 'N2', 'CO2', 'SO2', 'H2S',
                             bst.Chemical('Ash', search_db=False, default=True, phase='s', MW=1)], cache=True)
    bst.F.set_flowsheet('boiler_without_lime')
    with bst.System('sys') as sys:
        H = bst.HXutility('H', bst.Stream('feed', Water=100), T=360)
        B = bst.Boiler('B', ins=[bst.Stream('solids', H2S=0.1)])
    with pytest.raises(RuntimeError):
        sys.simulate()

def test_refrigeration_package_not_ready():
    with pytest.raises(NotImplementedError):
        bst.facilities.RefrigerationPackage('RP')

def test_chemical_capital_investment():
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    bst.F.set_flowsheet('chemical_capital_investment')
    with bst.System('sys') as sys:
        T = bst.StorageTank('T', bst.Stream('feed', Water=100, Ethanol=10), tau=1)
        CCI = bst.facilities.ChemicalCapitalInvestment('CCI', 'Ethanol', 1.)
    sys.simulate()
    assert CCI.design_results['Volume'] > 0
    assert CCI.purchase_cost > 0

def test_process_water_center_prices():
    bst.settings.set_thermo(['Water'], cache=True)
    bst.F.set_flowsheet('process_water_center_coverage')
    RO_price = bst.stream_utility_prices['Reverse osmosis water']
    process_price = bst.stream_utility_prices['Process water']
    try:
        makeup = bst.Stream('makeup', Water=10)
        process = bst.Stream('process', Water=10)
        PWC = bst.ProcessWaterCenter(
            'PWC', makeup_water_streams=[makeup], 
            process_water_streams=[makeup, process],
            reverse_osmosis_water_price=0.001, process_water_price=0.0005,
        )
        assert PWC.process_water_streams == [process]
        assert PWC.reverse_osmosis_water_price == 0.001
        assert PWC.process_water_price == 0.0005
        assert PWC.recycled_reverse_osmosis_grade_water is PWC.ins[0]
        assert PWC.recycled_process_water is PWC.ins[2]
        assert PWC.makeup_process_water is PWC.ins[3]
    finally:
        bst.stream_utility_prices['Reverse osmosis water'] = RO_price
        bst.stream_utility_prices['Process water'] = process_price

if __name__ == '__main__':
    test_boiler()
    test_boiler_turbogenerator_simple()
    test_boiler_without_lime()
    test_refrigeration_package_not_ready()
    test_chemical_capital_investment()
    test_process_water_center_prices()
