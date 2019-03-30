# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:48:21 2018

@author: yoelr
"""
from biosteam import Unit
from biosteam.units.splitter import Splitter


class Separator(Unit):
    _N_outs = 2
    _kwargs = Splitter._kwargs
    _run = Splitter._run
