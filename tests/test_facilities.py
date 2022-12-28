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
from biorefineries import cane

def test_facility_inheritance():
    with pytest.raises(bst.exceptions.UnitInheritanceError):
        class NewFacility(bst.Facility): pass
    
    class NewFacility(bst.Facility): 
        network_priority = 2


def test_boiler_turbogenerator():
    chemicals = cane.create_sugarcane_chemicals()
    chemicals.define_group(
        name='Fiber',
        IDs=['Cellulose', 'Hemicellulose', 'Lignin'],
        composition=[0.4704 , 0.2775, 0.2520],
        wt=True, # Composition is given as weight
    )
    bst.settings.set_thermo(chemicals)
    dilute_ethanol = bst.Stream('dilute_ethanol', Water=1390, Ethanol=590)
    bagasse = bst.Stream('bagasse', Water=0.4, Fiber=0.6, total_flow=8e4, units='kg/hr')
    with bst.System('sys') as sys:
        D1 = bst.BinaryDistillation('D1', ins=dilute_ethanol, Lr=0.999, Hr=0.89, k=1.25, LHK=('Ethanol', 'Water'))
        BT = bst.BoilerTurbogenerator('BT')
        BT.ins[0] = bagasse
    sys.simulate()
    assert_allclose(
        -BT.results().loc['Low pressure steam', 'Duty']['BT'],
        D1.results().loc['Low pressure steam', 'Duty']['D1'],
    )
    assert sys.power_utility.rate < 0. # Make sure there is excess electricity

    
if __name__ == '__main__':
    test_facility_inheritance()
    test_boiler_turbogenerator()