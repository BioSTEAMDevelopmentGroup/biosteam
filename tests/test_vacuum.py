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

def test_rotary_vane_vacuum_pump():

    bst.settings.set_thermo(['Water', 'CH4'])
    gas = bst.Stream('gas', CH4=2, Water=0.1)
    vac = bst.VacuumSystem(F_mass=gas.F_mass, F_vol=gas.F_vol, 
        P_suction=31325., vessel_volume=1.9e-3)
    assert "Rotary-vane pump, one stage" in vac.baseline_purchase_costs
    #!!! TODO: add test to purchase cost after updating CEPCI
    
    vac = bst.VacuumSystem(F_mass=gas.F_mass*10, F_vol=gas.F_vol*10, 
        P_suction=31325., vessel_volume=1.9e-3)
    assert "Rotary-vane pump, two stage" in vac.baseline_purchase_costs
    pass

if __name__ == '__main__':
    test_rotary_vane_vacuum_pump()


