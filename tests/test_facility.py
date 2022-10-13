# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose

def test_facility_inheritance():
    
    with pytest.raises(bst.exceptions.UnitInheritanceError):
        class NewFacility(bst.Facility): pass
    
    class NewFacility(bst.Facility): 
        network_priority = 2
    
if __name__ == '__main__':
    test_facility_inheritance()