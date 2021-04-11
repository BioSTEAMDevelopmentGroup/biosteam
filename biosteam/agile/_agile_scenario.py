# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains the abstract AgileScenario class that serves to create
objects that can retrive general results from multiple system scenarios
in such a way that it represents an agile production process.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.agile.AgileScenario

"""

__all__ = ('AgileScenario',)

class AgileScenario:
    """
    Abstract AgileScenario class for creating objects which may serve to retrive
    general results from multiple system scenarios in such a way that it 
    represents an agile production process. AgileScenario subclasses must be able 
    to create scenarios from systems and compile them to retrieve results later.
    
    """
    def __init_subclass__(cls):
        if not hasattr(cls, 'compile_scenarios'):
            raise NotImplementedError("missing method 'compile_scenarios'")
        if not hasattr(cls, 'create_scenario'):
            raise NotImplementedError("missing method 'create_scenario'")    
