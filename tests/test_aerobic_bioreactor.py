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
from biorefineries.sugarcane import chemicals
from numpy.testing import assert_allclose
import pytest

def test_jacketed_bubble_column():
    # Get arguments ready
    bst.settings.set_thermo(chemicals)
    feed = bst.Stream(
        Water=1.20e+05,
        Glucose=2.5e+04,
        units='kg/hr',
        T=32+273.15
    )
    rxn = bst.Rxn(
        'Glucose + O2 -> H2O + CO2', reactant='Glucose', X=0.5,
        correct_atomic_balance=True
    ) 
    kwargs = dict(
        ins=[feed, bst.Stream('air', phase='g')],
        outs=('vent', 'product'), tau=12, V_max=500,
        design='Bubble column', reactions=rxn,
    )
    
    # Must raise error if optimized power because default 
    # kLa method does not depend on agitation.
    R1 = bst.AeratedBioreactor(optimize_power=True, **kwargs)
    with pytest.raises(NotImplementedError):
        R1.simulate()

    R1 = bst.AeratedBioreactor(**kwargs)
    assert not R1.optimize_power, 'optimize power should default to false'
    
    # Test design results
    R1.simulate()
    expected_results = {
        'Reactor volume': 460.9817488148036,
        'Batch time': 18.30995695667444,
        'Loading time': 3.3099569566744393,
        'Number of reactors': 6,
        'Residence time': 12,
        'Length': 57.138719611650345,
        'Diameter': 19.046239870550114,
        'Weight': 74806.85,
        'Wall thickness': 0.39430441043892456
    }
    for name, expected in expected_results.items():
        assert_allclose(R1.design_results[name], expected)

    assert R1.design_results['Vessel type'] == 'Vertical', 'vessel type must be vertial'

    # Test gas flow rate sensitivity to aspect ratio
    air_flow_rates = []
    for L2D in (1, 3, 5):
        R1.length_to_diameter = L2D
        R1.simulate()
        air_flow_rates.append(R1.air.F_mass)
    
    assert_allclose(
        air_flow_rates, 
        [387392.31060733745, 
         191184.84862814427, 
         139537.07935482083]
    )
    
    # Must raise error if titer/yield is not feasible (mass transfer limitation)
    R1.length_to_diameter = 0.01
    with pytest.raises(RuntimeError): R1.simulate()


def test_jacketed_stirred_tank():
    # Get arguments ready
    bst.settings.set_thermo(chemicals)
    feed = bst.Stream(
        Water=1.20e+05,
        Glucose=2.5e+04,
        units='kg/hr',
        T=32+273.15
    )
    rxn = bst.Rxn(
        'Glucose + O2 -> H2O + CO2', reactant='Glucose', X=0.5,
        correct_atomic_balance=True
    ) 
    kwargs = dict(
        ins=[feed, bst.Stream('air', phase='g')],
        outs=('vent', 'product'), tau=12, V_max=500,
        reactions=rxn,
    )
    R1 = bst.AeratedBioreactor(**kwargs)
    assert R1.design == 'Stirred tank'
    assert R1.optimize_power, 'optimize power should default to True'
    
    # Test design results
    R1.simulate()
    expected_results = {
        'Reactor volume': 460.9817488148036,
        'Batch time': 18.30995695667444,
        'Loading time': 3.3099569566744393,
        'Number of reactors': 6,
        'Residence time': 12,
        'Length': 57.138719611650345,
        'Diameter': 19.046239870550114,
        'Weight': 74806.85,
        'Wall thickness': 0.39430441043892456
    }
    for name, expected in expected_results.items():
        assert_allclose(R1.design_results[name], expected)

    assert R1.design_results['Vessel type'] == 'Vertical', 'vessel type must be vertial'

    # Test gas flow rate sensitivity to aspect ratio
    air_flow_rates = []
    for L2D in (1, 3, 5):
        R1.length_to_diameter = L2D
        R1.simulate()
        air_flow_rates.append(R1.air.F_mass)
    
    assert_allclose(
        air_flow_rates, 
        [220615.2844394027,
         140308.29890959567, 
         116120.50925172029]
    )
    

if __name__ == '__main__':
    test_jacketed_bubble_column()
    test_jacketed_stirred_tank()
    