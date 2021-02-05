# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam.utils import format_title
__all__ = ('element_name',)

def element_name(element):
    if element:
        if isinstance(element, str):
            return element.replace('_', ' ')
        elif hasattr(element, 'line'):
            return element.line + '-' + element.ID.replace('_', ' ')
        else:
            unformatted_name = (element if isinstance(element, type) else type(element)).__name__
            return format_title(unformatted_name)
    else:
        return ''