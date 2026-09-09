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
import numpy as np
from numpy.testing import assert_allclose
import pytest

def test_gas_fed_bioreactor():
    bst.settings.set_thermo([
        'H2', 'CO2', 'CO', 'N2', 'O2', 'H2O', 'AceticAcid'
    ])

    # Feeds
    media = bst.Stream(H2O=20e3, units='kg/hr')
    H2 = bst.Stream(H2=1e3, units='kg/hr', phase='g')
    fluegas = bst.Stream(CO=23, CO2=23, H2=4.5, N2=49.5, units='kg/hr', phase='g') 
    vent = bst.Stream(phase='g')
    effluent = bst.Stream(phase='l')

    # Model acetic acid production from H2 and CO2
    substrate_reaction = bst.Reaction(
        'CO + H2O -> CO2 + H2',
        reactant='CO', correct_atomic_balance=True, X=1
    )
    rxn = bst.Reaction(
        'H2 + CO2 -> AceticAcid + H2O',
        reactant='CO2', correct_atomic_balance=True, X=1
    )
    
    def assert_titer():
        titer = sum([i.imass['AceticAcid'] for i in bioreactor.outs]) / media.F_mass * 1000
        assert_allclose(titer, bioreactor.titer['AceticAcid'], rtol=1e-3)
    
    def assert_flow_rates_the_same():
        assert_allclose(flow_rates_baseline, get_flow_rates(), rtol=1e-3)
    
    def get_flow_rates():
        return np.array([media.F_mass, H2.F_mass, fluegas.F_mass])
    
    ### Solve H2 and Flue gas flow rates (2 gas stream) ###
    H2.F_mass = 1
    fluegas.F_mass = 1
    bioreactor = bst.GasFedBioreactor(
        ins=[media, H2, fluegas], 
        outs=(vent, effluent), 
        controlled_feeds=[2, 1],
        tau=68, V_max=500,
        length_to_diameter=12,
        reactions=rxn, 
        substrate_reactions=substrate_reaction,
        gas_substrates=('H2', 'CO', 'CO2'),
        controlled_gas_substrates=('CO2', 'H2'),
        design='Bubble column',
        kW_per_m3=0,
        batch=False,
        T=305.15,
        titer={'AceticAcid': 20},
    )
    bioreactor.simulate()
    flow_rates_baseline = get_flow_rates()
    assert_titer()
    
    ### Solve media and H2 (1 liquid stream, 1 gas stream) ###
    media.F_mass = 1
    H2.F_mass = 1
    bioreactor = bst.GasFedBioreactor(
        ins=[media, H2, fluegas], 
        outs=(vent, effluent), 
        controlled_feeds=[0, 1],
        tau=68, V_max=500,
        length_to_diameter=12,
        reactions=rxn, 
        substrate_reactions=substrate_reaction,
        gas_substrates=('H2', 'CO', 'CO2'),
        controlled_gas_substrates=('CO2', 'H2'),
        design='Bubble column',
        kW_per_m3=0,
        batch=False,
        T=305.15,
        titer={'AceticAcid': 20},
    )
    bioreactor.simulate()
    assert_titer()
    assert_flow_rates_the_same()

    ### Solve media (1 liquid stream) ###
    media.F_mass = 1
    bioreactor = bst.GasFedBioreactor(
        ins=[media, H2, fluegas], 
        outs=(vent, effluent), 
        controlled_feeds=[0],
        tau=68, V_max=500,
        length_to_diameter=12,
        reactions=rxn, 
        substrate_reactions=substrate_reaction,
        gas_substrates=('H2', 'CO', 'CO2'),
        design='Bubble column',
        kW_per_m3=0,
        batch=False,
        T=305.15,
        titer={'AceticAcid': 20},
    )
    bioreactor.simulate()
    assert_titer()
    assert_flow_rates_the_same()
    

if __name__ == '__main__':
    test_gas_fed_bioreactor()
    