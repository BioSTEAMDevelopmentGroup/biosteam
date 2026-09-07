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
from numpy.testing import assert_allclose
import pytest

def test_gas_fed_bioreactor():
    bst.settings.set_thermo(['H2', 'CO2', 'CO', 'N2', 'O2', 'H2O', 'AceticAcid'])

    # Feeds
    media = bst.Stream(ID='media', H2O=20e3, units='kg/hr')
    H2 = bst.Stream(ID='H2', H2=100, units='kg/hr', phase='g')
    fluegas = bst.Stream(CO=23, CO2=23, H2=4.5, N2=49.5, units='m3/hr', phase='g') 

    # Model acetic acid production from H2 and CO2
    # Model acetic acid production from H2 and CO2
    substrate_reaction = bst.Reaction(
        'CO + H2O -> CO2 + H2',
        reactant='CO', correct_atomic_balance=True, X=1
    )
    rxn = bst.Reaction(
        'H2 + CO2 -> AceticAcid + H2O',
        reactant='CO2', correct_atomic_balance=True, X=1
    )
    bioreactor = bst.GasFedBioreactor(
        ins=[media, H2, fluegas], 
        outs=('vent', 'product'), 
        tau=68, V_max=500,
        length_to_diameter=6,
        reactions=rxn, 
        substrate_reactions=substrate_reaction,
        gas_substrates=('H2', 'CO', 'CO2'),
        design='Bubble column',
        kW_per_m3=0,
        batch=False,
        T=305.15
    )
    bioreactor.simulate()
    

if __name__ == '__main__':
    test_gas_fed_bioreactor()
    