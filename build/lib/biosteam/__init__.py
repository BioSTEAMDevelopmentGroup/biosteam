#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:28:09 2017

@author: Yoel Rene Cortes-Pena
"""
__all__ = ['Species', 'WorkingSpecies', 'Stream', 'MixedStream',
            'Unit', 'System', 'TEA', 'PowerUtility', 'HeatUtility',
            'find', 'Flowsheet', 'CE']

from lazypkg import LazyPkg
LazyPkg(__name__, ['_equilibrium', '_utils', 'units',
                   'evaluation', 'inspect', 'compounds',
                   'reaction'])

#: Chemical engineering plant cost index (defaults to 567.5 at 2017)
CE = 567.5 

# %% Import base utils

import pandas as pd
import numpy as np
from pint import UnitRegistry
import os

# Set pint Unit Registry
_ureg = UnitRegistry()
_ureg.default_format = '~P'
_ureg.load_definitions(os.path.dirname(os.path.realpath(__file__)) + '/my_units_defs.txt')
_Q = _ureg.Quantity
_Q._repr_latex_ = _Q._repr_html_ = _Q.__str__ = _Q.__repr__ = lambda self: self.__format__('')

# Set number of digits displayed
np.set_printoptions(suppress=False)
np.set_printoptions(precision=3) 
pd.options.display.float_format = '{:.3g}'.format
pd.set_option('display.max_rows', 35)
pd.set_option('display.max_columns', 10)
pd.set_option('max_colwidth', 35)
del np, pd, os, UnitRegistry


# %% Initialize BioSTEAM 

from ._species import Species, WorkingSpecies
from ._stream import Stream
from ._mixed_stream import MixedStream
from ._heat_utility import HeatUtility
from ._power_utility import PowerUtility
from ._unit import Unit
from ._system import System
from ._tea import CombinedTEA, TEA
from ._flowsheet import Flowsheet, find




