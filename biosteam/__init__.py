#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:28:09 2017

@author: Yoel Rene Cortes-Pena
"""

__all__ = ['units', 'evaluation', 'inspect', 'find', 'compounds', 'Species', 'Stream', 'MixedStream', 'HeatUtility', 'PowerUtility', 'Unit', 'System', 'TEA']


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

# %% Import biosteam classes

from ._species import *
from ._stream import *
from ._mixed_stream import *
from ._heat_utility import *
from ._power_utility import *
from ._unit import *
from ._system import *
from ._tea import *
from ._flowsheet import *

from . import inspect
from . import compounds
from . import _species 
from . import _stream
from . import _mixed_stream
from . import _unit
from . import _heat_utility
from . import _power_utility
from . import evaluation
from . import _system 
from . import _tea
from . import _flowsheet
from . import units

__all__.extend(_species.__all__)
__all__.extend(_stream.__all__)
__all__.extend(_mixed_stream.__all__)
__all__.extend(_heat_utility.__all__)
__all__.extend(_power_utility.__all__)
__all__.extend(_unit.__all__)
__all__.extend(_system.__all__)
__all__.extend(_tea.__all__)
__all__.extend(_flowsheet.__all__)
