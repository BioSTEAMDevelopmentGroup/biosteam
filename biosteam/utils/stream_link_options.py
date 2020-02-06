# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 00:45:47 2019

@author: yoelr
"""
from collections import namedtuple
from .not_implemented_method import NotImplementedMethod

__all__ = ('StreamLinkOptions',
           'static',
           'static_flow',
           'static_flow_and_phase',
           'static_link_options',
           'static_flow_link_options',
           'static_flow_and_phase_options')

# %% Linking options

StreamLinkOptions = namedtuple('StreamLinkOptions', ('flow', 'phase', 'TP'), module=__name__)
static_link_options = StreamLinkOptions(flow=True, TP=True, phase=True)
static_flow_and_phase_options = StreamLinkOptions(flow=True, TP=False, phase=True)
static_flow_link_options = StreamLinkOptions(flow=True, TP=False, phase=False)

def static(cls):
    cls._stream_link_options = static_link_options
    cls._run = NotImplementedMethod
    cls._N_ins = cls._N_outs = 1
    return cls

def static_flow(cls):
    cls._stream_link_options = static_flow_link_options
    cls._N_ins = cls._N_outs = 1
    return cls

def static_flow_and_phase(cls):
    cls._stream_link_options = static_flow_and_phase_options
    cls._N_ins = cls._N_outs = 1
    return cls

del namedtuple