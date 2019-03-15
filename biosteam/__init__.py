#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:28:09 2017

@author: Yoel Rene Cortes-Pena
"""
name = 'biosteam'

__all__ = ['exceptions', 'utils', 'find', 'Compound', 'Chemical', 'DissolvedCompound', 'Species', 'Stream', 'MixedStream', 'HeatUtility', 'PowerUtility', 'Unit', 'System', 'ureg', 'Quantity', 'Q_', 'units_of_measure', 'TEA', 'ProxyStream']


# %% Import base utils

import pandas as pd
import numpy as np
from bookkeep.unit_registry import ureg, Q_
from bookkeep import ReadOnlyBook, SmartBook

# Allias
Quantity = Q_

# Set number of digits displayed
np.set_printoptions(suppress=False)
np.set_printoptions(precision=3) 
pd.options.display.float_format = '{:.3g}'.format
pd.set_option('display.max_rows', 35)
pd.set_option('display.max_columns', 10)
pd.set_option('max_colwidth', 35)

# Biosteam units of measure
units_of_measure = ReadOnlyBook(MW='g/mol',
                                mass='kg/hr',
                                mol='kmol/hr',
                                vol='m^3/hr',
                                massnet='kg/hr',
                                molnet='kmol/hr',
                                volnet='m^3/hr',
                                massfrac='kg/kg',
                                molfrac='kmol/kmol',
                                volfrac='m^3/m^3',
                                T='K',
                                P='Pa',
                                H='kJ/hr',
                                S='kJ/hr',
                                G='kJ/hr',
                                U='kJ/hr',
                                A='kJ/hr',
                                Hf='kJ/hr',
                                C='kJ/K/hr',
                                Vm='m^3/mol',
                                Cpm='J/mol/K',
                                Cp='J/g/K',
                                rho='kg/m^3',
                                rhom='mol/m^3',
                                nu='m^2/s',
                                mu='Pa*s',
                                sigma='N/m',
                                k='W/m/K',
                                alpha='m^2/s')


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
from .reports import *
from .plots import *
from .units import *

from . import exceptions
from . import utils
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
from . import reports
from . import plots
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
__all__.extend(reports.__all__)
__all__.extend(plots.__all__)
__all__.extend(units.__all__)
SmartBook.Warning = exceptions.DesignWarning
