# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autodata:: stream_utility_prices
.. autodata:: CE

"""
from __future__ import annotations
__version__ = '2.33.3'

#: Chemical engineering plant cost index (defaults to 567.5 at 2017).
CE: float = 567.5 

#: User-defined impact indicators and their units of measure.
impact_indicators: dict[str, str] = {}

#: Price of stream utilities [USD/kg] which are defined as 
#: inlets and outlets to unit operations.
stream_utility_prices: dict[str, float] = {'Natural gas': 0.218,
                                           'Ash disposal': -0.0318}

# %% Workaround for readthedocs, which fails to cache numba

import numba
try:
    @numba.njit(cache=True)
    def f_dummy(): pass
except RuntimeError: # pragma: no cover
    def njit(*args, **kwargs):
        kwargs['cache'] = False
        return numba.jit(*args, **kwargs)
    numba.njit = njit
    del njit
else:
    del f_dummy

# %% Initialize BioSTEAM 

from thermosteam import Chemical, Chemicals, Thermo, Stream, MultiStream, settings, ProcessSettings, speed_up
from ._preferences import preferences
from ._heat_utility import *
from ._power_utility import PowerUtility
from . import plots
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
from . import hidden_connection
from . import units
from .units import *
from . import evaluation
from .evaluation import *
from . import exceptions
from . import report

__all__ = (
    'Unit', 'PowerUtility', 'UtilityAgent', 'HeatUtility',
    'utils', 'units', 'evaluation', 'Chemical', 'Chemicals', 'Stream',
    'MultiStream', 'settings', 'exceptions', 'report',
    'process_tools', 'preferences', *_system.__all__, *_flowsheet.__all__, 
    *_tea.__all__, *units.__all__, *evaluation.__all__, 
    *process_tools.__all__, *hidden_connection.__all__,
)

from .hidden_connection import *

def nbtutorial():
    preferences.reset()
    preferences.light_mode(bg='#ffffffaa')
    preferences.tooltips_full_results = False
    preferences.graphviz_format = 'html'
    from warnings import filterwarnings
    filterwarnings('ignore')