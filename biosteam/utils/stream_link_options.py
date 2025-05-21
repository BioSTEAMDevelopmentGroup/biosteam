# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from collections import namedtuple

__all__ = ('static',
           'static_flow_and_phase')

# %% Linking options


def _run_static(self):
    self._outs[0].copy_like(self._ins[0])

def static(cls=None):
    if cls is None: return lambda cls: static(cls)
    cls._link_streams = True
    cls._run = _run_static
    return cls

def _run_static_flow(self):
    outlet = self._outs[0]
    inlet = self._ins[0]
    outlet.phases = inlet.phases
    outlet.copy_flow(inlet)

def static_flow_and_phase(cls):
    cls._N_ins = cls._N_outs = 1
    if '_run' not in cls.__dict__: cls._run = cls._run_static_flow
    return cls

del namedtuple