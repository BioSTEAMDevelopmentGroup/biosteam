# -*- coding: utf-8 -*-
<<<<<<< HEAD
"""
Created on Tue May 14 14:35:55 2019

@author: yoelr
=======
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
>>>>>>> cd2c5013aaf9b5bc94bb764b52fd37db183472f1
"""
import re
__all__ = ('element_name',)

def element_name(element):
    if element:
        if isinstance(element, type):
            return re.sub(r"\B([A-Z])", r" \1", element.__name__.replace('_', ' ')).capitalize()
        elif isinstance(element, str):
            return element.replace('_', ' ')
        else:
            return element.line + '-' + element.ID.replace('_', ' ')
    else:
        return ''