# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest

def test_model():
    from biosteam.examples import ethanol_subsystem as ethanol_sys
    from biosteam import main_flowsheet as F
    import biosteam as bst
    
    # Make sure repeated metrics raise an error
    biorefinery = bst.process_tools.UnitGroup('Biorefinery', ethanol_sys.units)
    heating_duty = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    heating_duty_repeated = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    with pytest.raises(ValueError):
        model = bst.Model(ethanol_sys, [heating_duty, heating_duty_repeated])
    
    model = bst.Model(ethanol_sys, [heating_duty])
    
    # Make sure repeated parameters raise an error
    R301 = F.flowsheet.ethanol_subsystem_example.unit.R301
    @model.parameter(element=R301)
    def set_efficiency(efficiency):
        R301.efficiency = efficiency
    
    with pytest.raises(ValueError):
        @model.parameter(element=R301)
        def set_efficiency(efficiency):
            R301.efficiency = efficiency
    
    bst.process_tools.default_utilities()
    bst.CE = 567.5