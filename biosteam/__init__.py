# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
__version__ = '2.38.3'

#: Chemical engineering plant cost index (defaults to 567.5 at 2017).
CE: float = 567.5 

#: User-defined impact indicators and their units of measure.
impact_indicators: dict[str, str] = {}

#: Price of stream utilities [USD/kg] which are defined as 
#: inlets and outlets to unit operations.
stream_utility_prices: dict[str, float] = {}

#: Defined allocation property and basis pairs for LCA.
allocation_properties: dict[str, str] = {}

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

import thermosteam
from thermosteam import *
from ._preferences import preferences
from ._heat_utility import UtilityAgent, HeatUtility
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
from . import _temporary_connection
from ._temporary_connection import *
from . import units
from .units import *
from ._facility import Facility
from . import facilities
from .facilities import *
from . import wastewater
from .wastewater import *
from . import evaluation
from .evaluation import *
from . import exceptions
from . import report
from . import _settings

__all__ = (
    'Unit', 'PowerUtility', 'UtilityAgent', 'HeatUtility', 'Facility',
    'utils', 'units', 'facilities', 'wastewater', 'evaluation', 'Chemical', 'Chemicals', 'Stream',
    'MultiStream', 'settings', 'exceptions', 'report',
    'process_tools', 'preferences', *_system.__all__, *_flowsheet.__all__, 
    *_tea.__all__, *units.__all__, *facilities.__all__, *wastewater.__all__,
    *evaluation.__all__, *process_tools.__all__, *_temporary_connection.__all__,
)

def nbtutorial():
    preferences.reset()
    preferences.light_mode(bg='#ffffffaa')
    preferences.tooltips_full_results = False
    preferences.graphviz_format = 'html'
    from warnings import filterwarnings
    filterwarnings('ignore')
