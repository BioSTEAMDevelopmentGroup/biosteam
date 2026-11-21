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

def test_utility_agent_properties():
    bst.settings.set_thermo(['Water'], cache=True)
    agent = bst.HeatUtility.get_agent('low_pressure_steam')
    assert agent.isheating_agent
    with pytest.raises(AttributeError): agent.price
    with pytest.raises(AttributeError): agent.price = 1
    with pytest.raises(AttributeError): agent.cost
    with pytest.raises(AttributeError): agent.cost = 1
    # Representation with all price combinations
    original = (agent.heat_transfer_price, agent.regeneration_price)
    try:
        for ht, rg in [(1e-6, 0), (0, 1e-3), (1e-6, 1e-3), (0, 0)]:
            agent.heat_transfer_price = ht
            agent.regeneration_price = rg
            assert 'price' in agent._info()
    finally:
        agent.heat_transfer_price, agent.regeneration_price = original
    assert 'Water' not in agent._info(N=0)
    df = agent._info(df=True)
    assert 'Temperature [K]' in df.index
    # Agent with T_limit and many chemicals
    bst.settings.set_thermo(['Water', 'Ethanol', 'Methanol'], cache=True)
    stream = bst.Stream('brine', Water=1, Ethanol=1, Methanol=1, T=250)
    new_agent = stream.to_agent(T_limit=260)
    assert new_agent.iscooling_agent
    assert '...' in new_agent._info(N=1)
    df = new_agent._info(N=1, df=True)
    assert 'Remainder [%]' in df.index
    empty = bst.UtilityAgent('empty_agent', Water=1, units='kg/hr')
    assert_allclose(empty.F_mol, 1)

def test_heat_utility_agent_lookup():
    bst.settings.set_thermo(['Water'], cache=True)
    hu = bst.HeatUtility()
    assert repr(hu) == '<HeatUtility: None>'
    with pytest.raises(LookupError): bst.HeatUtility.get_agent('none')
    with pytest.raises(LookupError): bst.HeatUtility.get_heating_agent('none')
    with pytest.raises(LookupError): bst.HeatUtility.get_cooling_agent('none')
    assert bst.HeatUtility.get_heating_agent('low_pressure_steam').ID == 'low_pressure_steam'
    assert bst.HeatUtility.get_cooling_agent('cooling_water').ID == 'cooling_water'
    with pytest.raises(RuntimeError): bst.HeatUtility.get_suitable_cooling_agent(1.)
    heating_agents = bst.HeatUtility.heating_agents
    try:
        bst.HeatUtility.heating_agents = [i for i in heating_agents if not i.isfuel]
        with pytest.raises(RuntimeError): bst.HeatUtility.get_suitable_heating_agent(1e4)
    finally:
        bst.HeatUtility.heating_agents = heating_agents
    with pytest.raises(ValueError):
        bst.HeatUtility.get_inlet_and_outlet_pinch_temperature(False, 400, 300)
    with pytest.raises(AttributeError): hu.inlet_process_stream
    with pytest.raises(AttributeError): hu.outlet_process_stream
    hu.ID = 'cooling_water'
    assert hu.agent.ID == 'cooling_water'

def test_heat_utility_set_by_flow_rate_and_scale():
    bst.settings.set_thermo(['Water'], cache=True)
    steam = bst.HeatUtility.get_agent('low_pressure_steam')
    hu = bst.HeatUtility()
    hu.set_utility_by_flow_rate(steam, 10)
    assert hu.duty > 0 and hu.cost > 0
    hu.set_utility_by_flow_rate(steam, 0)
    assert not hu
    # Load same agent twice and then a different one
    hu(1e6, 300, 350)
    hu(2e6, 300, 350)
    assert hu.duty > 0
    hu(-1e6, 350, 320)
    assert hu.duty < 0
    with pytest.raises(ValueError):
        hu.set_utility_by_flow_rate(bst.HeatUtility.get_agent('cooling_water'), 10)
    # Fuel agent
    fuel = [i for i in bst.HeatUtility.heating_agents if i.isfuel]
    if fuel:
        hu = bst.HeatUtility()
        hu(1e6, 300, 700)
        assert hu.agent.isfuel
        # Note: hu.scale fails for fuel agents (read-only outlet flow)
        hu.empty()
        hu2 = bst.HeatUtility()
        hu2(1e6, 300, 700) # Reuses dumped fuel streams
        assert hu2.duty > 0
    # Mix from with a single utility
    hu = bst.HeatUtility()
    other = bst.HeatUtility()
    other(1e6, 300, 350)
    hu.mix_from([other])
    assert_allclose(hu.duty, other.duty)

def test_heat_utility_characterization_factors():
    bst.settings.set_thermo(['Water'], cache=True)
    CFs = bst.HeatUtility.characterization_factors.copy()
    indicators = bst.settings.impact_indicators.copy()
    try:
        bst.settings.define_impact_indicator('GWP', 'kg*CO2e')
        bst.HeatUtility.set_CF('low_pressure_steam', 'GWP', 1., basis='kmol')
        bst.HeatUtility.set_CF('cooling_water', 'GWP', -1e-6, basis='kJ', units='g*CO2e')
        bst.HeatUtility.set_CF('high_pressure_steam', 'GWP', 0.1, basis='kg')
        with pytest.raises(ValueError):
            bst.HeatUtility.set_CF('cooling_water', 'GWP', 1e-6, basis='kJ')
        with pytest.raises(Exception):
            bst.HeatUtility.set_CF('cooling_water', 'GWP', 1e-6, basis='m')
        assert bst.HeatUtility.get_CF('none', 'GWP', basis='kg') == 0.
        assert bst.HeatUtility.get_CF('low_pressure_steam', 'GWP', units='g*CO2e')[0] > 0
        # Note: get_impact/get_inventory fail for molar basis ('kmol' != 'mol')
        hu = bst.HeatUtility()
        hu(-1e6, 350, 320, agent=bst.HeatUtility.get_agent('cooling_water'))
        assert hu.get_impact('GWP') > 0
        assert hu.get_inventory('GWP')[1] == 'kJ/hr'
        assert hu.has_impact('GWP')
        hu(1e6, 300, 450, agent=bst.HeatUtility.get_agent('high_pressure_steam'))
        assert hu.get_impact('GWP') > 0
        assert hu.get_inventory('GWP')[1] == 'kg/hr'
        hu(1e6, 300, 350, agent=bst.HeatUtility.get_agent('medium_pressure_steam'))
        assert hu.get_impact('GWP') == 0.
    finally:
        bst.HeatUtility.characterization_factors = CFs
        bst.settings.impact_indicators = indicators

def test_power_utility_characterization_factors():
    CFs = bst.PowerUtility.characterization_factors.copy()
    indicators = bst.settings.impact_indicators.copy()
    try:
        bst.settings.define_impact_indicator('GWP', 'kg*CO2e')
        bst.PowerUtility.set_CF('GWP', 0.5, 0.4, units='g*CO2e')
        assert bst.PowerUtility.get_CF('GWP', production=False) > 0
        assert bst.PowerUtility.get_CF('GWP', consumption=False) > 0
        assert bst.PowerUtility.get_CF('GWP', False, False) is None
        assert bst.PowerUtility.get_CF('GWP', units='g*CO2e', production=False) > 0
        assert len(bst.PowerUtility.get_CF('GWP', units='g*CO2e')) == 2
        assert bst.PowerUtility.get_CF('GWP', basis='MWhr', production=False) is not None
        bst.PowerUtility.set_CF('GWP', consumption=0.5)
        bst.PowerUtility.set_CF('GWP', production=0.5)
        with pytest.raises(ValueError):
            bst.PowerUtility.set_CF('GWP')
        pu = bst.PowerUtility(consumption=10)
        assert pu.get_impact('GWP') > 0
        assert pu.has_impact('GWP')
        assert bst.PowerUtility().get_impact('none') == 0.
        assert not bst.PowerUtility().has_impact('GWP')
        pu2 = pu.copy()
        pu2(5)
        assert pu2.power == 5
    finally:
        bst.PowerUtility.characterization_factors = CFs
        bst.settings.impact_indicators = indicators

def test_settings_properties():
    settings = bst.settings
    CE = settings.CEPCI
    cooling = settings.cooling_agents
    heating = settings.heating_agents
    prices = settings.stream_prices
    indicators = settings.impact_indicators
    CFs = settings.utility_characterization_factors
    skip = settings.skip_simulation_of_units_with_empty_inlets
    try:
        settings.CEPCI = 600.
        assert bst.CE == 600.
        settings.cooling_agents = cooling
        settings.heating_agents = heating
        settings.stream_prices = prices
        settings.impact_indicators = indicators
        settings.utility_characterization_factors = CFs
        settings.skip_simulation_of_units_with_empty_inlets = True
        assert settings.skip_simulation_of_units_with_empty_inlets
        assert isinstance(settings.allocation_properties, dict)
        with pytest.raises(ValueError):
            settings.get_impact_indicator_units('not an indicator')
    finally:
        settings.CEPCI = CE
        settings.skip_simulation_of_units_with_empty_inlets = skip

if __name__ == '__main__':
    test_utility_agent_properties()
    test_heat_utility_agent_lookup()
    test_heat_utility_set_by_flow_rate_and_scale()
    test_heat_utility_characterization_factors()
    test_power_utility_characterization_factors()
    test_settings_properties()
