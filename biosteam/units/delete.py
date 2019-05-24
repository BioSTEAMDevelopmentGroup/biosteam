# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:30:37 2019

@author: yoelr
"""

from biosteam import units, Unit

splitter = units.metaclasses.splitter
class SubUnit(Unit, metaclass=splitter):
    pass