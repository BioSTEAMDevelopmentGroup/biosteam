# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from collections import namedtuple

__all__ = ('StreamLinkOptions',
           'static',
           'static_flow_and_phase',
           'static_link_options',
           'static_flow_link_options',
           'static_flow_and_phase_options')

# %% Linking options

StreamLinkOptions = namedtuple('StreamLinkOptions', ('flow', 'phase', 'TP'), module=__name__)
static_link_options = StreamLinkOptions(flow=True, TP=True, phase=True)
static_flow_and_phase_options = StreamLinkOptions(flow=True, TP=False, phase=True)
static_flow_link_options = StreamLinkOptions(flow=True, TP=False, phase=False)

def _run_static(self):
    self._outs[0].copy_like(self._ins[0])

def static(cls):
    cls._stream_link_options = static_link_options
    cls._N_ins = cls._N_outs = 1
    cls._run = _run_static
    return cls

def _run_static_flow(self):
    outlet = self._outs[0]
    inlet = self._ins[0]
    outlet.phases = inlet.phases
    outlet.copy_flow(inlet)

def static_flow_and_phase(cls):
    cls._stream_link_options = static_flow_and_phase_options
    cls._N_ins = cls._N_outs = 1
    return cls

del namedtuple