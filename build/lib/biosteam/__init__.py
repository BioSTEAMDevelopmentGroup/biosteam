#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:28:09 2017

@author: Yoel Rene Cortes-Pena
"""
name = 'biosteam'

__all__ = ['exceptions', 'utils', 'inspect', 'find', 'Compound', 'Chemical', 'DissolvedCompound', 'Species', 'Stream', 'MixedStream', 'HeatUtility', 'PowerUtility', 'Unit', 'System', 'TEA', 'ProxyStream', 'inspect', 'Sensitivity']


# %% Import base utils

import pandas as pd
import numpy as np
from pint import UnitRegistry
import os

# Remove latex formating
def _new_format(self): return self.__format__('')

# Set pint Unit Registry
ureg = UnitRegistry()
ureg.default_format = '~P'
dir_path = os.path.dirname(os.path.realpath(__file__)) + '/'
ureg.load_definitions(dir_path + 'my_units_defs.txt')
Q_ = ureg.Quantity
Q_._repr_latex_ = Q_._repr_html_ = Q_.__str__ = Q_.__repr__ = _new_format

# Set number of digits displayed
np.set_printoptions(suppress=False)
np.set_printoptions(precision=3) 
pd.options.display.float_format = '{:.3g}'.format
pd.set_option('display.max_rows', 35)
pd.set_option('display.max_columns', 10)
pd.set_option('max_colwidth', 35)


# %% Import biosteam classes

from .compound import *
from .chemical import *
from .dissolved_compound import *
from .species import *
from .stream import *
from .mixed_stream import *
from .proxy_stream import *
from .heat_utility import *
from .power_utility import *
from .unit import *
from .sim import *
from .system import *
from .tea import *
from .flowsheet import *
from .report import *
from .units import *

from . import exceptions
from . import utils
from . import inspect
from . import compound 
from . import chemical
from . import dissolved_compound 
from . import species 
from . import stream
from . import mixed_stream
from . import proxy_stream
from . import unit
from . import heat_utility
from . import power_utility
from . import sim
from . import system 
from . import tea
from . import flowsheet
from . import report
from . import units

__all__.extend(compound.__all__)
__all__.extend(chemical.__all__)
__all__.extend(dissolved_compound.__all__)
__all__.extend(species.__all__)
__all__.extend(stream.__all__)
__all__.extend(mixed_stream.__all__)
__all__.extend(proxy_stream.__all__)
__all__.extend(heat_utility.__all__)
__all__.extend(power_utility.__all__)
__all__.extend(unit.__all__)
__all__.extend(sim.__all__)
__all__.extend(system.__all__)
__all__.extend(tea.__all__)
__all__.extend(flowsheet.__all__)
__all__.extend(report.__all__)
__all__.extend(units.__all__)
