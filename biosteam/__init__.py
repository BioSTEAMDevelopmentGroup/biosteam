# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
__version__ = '2.32.7'

#: dict[str, float] Price of stream utilities (in USD/kg) which are defined as 
#: inlets and outlets to unit operations.
stream_utility_prices = {'Natural gas': 0.218,
                         'Ash disposal': -0.0318}

#: Chemical engineering plant cost index (defaults to 567.5 at 2017).
CE = 567.5 

#: Whether to label the ID of streams with sources and sinks in process 
#: flow diagrams.
LABEL_PROCESS_STREAMS_IN_DIAGRAMS = True

#: Whether to automatically generate diagrams when displaying an object in the
#: IPython console.
ALWAYS_DISPLAY_DIAGRAMS = True

#: Whether to ignore unit graphics and display unit nodes as dots in process
#: flow diagrams.
MINIMAL_UNIT_DIAGRAMS = False

#: Whether to number unit operations in diagrams according to their order in the system path.
LABEL_PATH_NUMBER_IN_DIAGRAMS = False

#: Whether to clock the simulation time of unit operations in diagrams.
PROFILE_UNITS_IN_DIAGRAMS = False

#: Whether to raise exception regarding problems displaying graphviz diagrams
RAISE_GRAPHVIZ_EXCEPTION = False

#: Graphviz background color
GRAPHVIZ_BACKGROUND_COLOR = 'transparent'

#: Color of streams in BioSTEAM graphviz diagrams
GRAPHVIZ_STREAM_COLOR = '#90918e'

#: Color of stream labels in BioSTEAM graphviz diagrams
GRAPHVIZ_LABEL_COLOR = '#90918e'

#: Color of subsystem clusters in BioSTEAM graphviz diagrams
GRAPHVIZ_DEPTH_COLORS = ('#7ac0836f',)

def _set_mode(stream, label, bg, cluster, save):
    global GRAPHVIZ_BACKGROUND_COLOR, GRAPHVIZ_STREAM_COLOR, GRAPHVIZ_LABEL_COLOR, GRAPHVIZ_DEPTH_COLORS
    GRAPHVIZ_BACKGROUND_COLOR = bg
    GRAPHVIZ_STREAM_COLOR = stream
    GRAPHVIZ_LABEL_COLOR = label
    GRAPHVIZ_DEPTH_COLORS = cluster
    if save: save_mode()

def classic_mode(stream=GRAPHVIZ_STREAM_COLOR, 
                 label=GRAPHVIZ_LABEL_COLOR, 
                 bg=GRAPHVIZ_BACKGROUND_COLOR,
                 cluster=GRAPHVIZ_DEPTH_COLORS,
                 save=True):
    _set_mode(stream, label, bg, cluster, save)

def dark_mode(stream='#98a2ad', label='#e5e5e5', bg='#000000bf', cluster=('#f3c3546f',), save=True):
    _set_mode(stream, label, bg, cluster, save)

def light_mode(stream='#4e4e4e', label='#90918e', bg='#ffffffdf', cluster=('#7ac083af',), save=True):
    _set_mode(stream, label, bg, cluster, save)

night_mode = dark_mode
day_mode = light_mode

def load_mode(file=None):
    if file is None:
        import os
        folder = os.path.dirname(__file__)
        file = os.path.join(folder, 'graphviz_color_settings.txt')
    try:
        global GRAPHVIZ_BACKGROUND_COLOR, GRAPHVIZ_STREAM_COLOR, GRAPHVIZ_LABEL_COLOR, GRAPHVIZ_DEPTH_COLORS
        with open(file) as f:
            (GRAPHVIZ_BACKGROUND_COLOR,
             GRAPHVIZ_STREAM_COLOR,
             GRAPHVIZ_LABEL_COLOR,
             GRAPHVIZ_DEPTH_COLORS) = f.readlines()
        GRAPHVIZ_DEPTH_COLORS = eval(GRAPHVIZ_DEPTH_COLORS)
    except:
        pass
    
def save_mode():
    import os
    folder = os.path.dirname(__file__)
    file = os.path.join(folder, 'graphviz_color_settings.txt')
    with open(file, 'wb') as f:
        f.write('\n'.join([GRAPHVIZ_BACKGROUND_COLOR, GRAPHVIZ_STREAM_COLOR, str(GRAPHVIZ_LABEL_COLOR)]).encode())

load_mode()

# %% Workaround for readthedocs, which fails to cache numba

import numba
try:
    @numba.njit(cache=True)
    def f_dummy(): pass
except RuntimeError:
    def njit(*args, **kwargs):
        kwargs['cache'] = False
        return numba.jit(*args, **kwargs)
    numba.njit = njit
    del njit
else:
    del f_dummy

# %% Initialize BioSTEAM 

from thermosteam import Chemical, Chemicals, Thermo, Stream, MultiStream, settings, speed_up
from ._heat_utility import HeatUtility, UtilityAgent
from ._power_utility import PowerUtility
from .utils import *
from ._unit import Unit
from . import _system
from ._system import *
from . import _flowsheet
from ._flowsheet import *
from . import process_tools
from .process_tools import *
from . import _tea
from ._tea import *
from . import utils
from . import units
from .units import *
from . import evaluation
from .evaluation import *
from . import exceptions
from . import report

__all__ = (
    'Unit', 'PowerUtility', 'HeatUtility', 'UtilityAgent',
    'utils', 'units', 'evaluation', 'Chemical', 'Chemicals', 'Stream',
    'MultiStream', 'settings', 'exceptions', 'report',
    'process_tools', *_system.__all__, *_flowsheet.__all__, 
    *_tea.__all__, *units.__all__, *evaluation.__all__, 
    *process_tools.__all__, 
)