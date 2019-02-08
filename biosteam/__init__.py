#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:28:09 2017

@author: Yoel Rene Cortes-Pena
"""
name = 'biosteam'

__all__ = ['exceptions', 'utils', 'find', 'Compound', 'Chemical', 'DissolvedCompound', 'Species', 'Stream', 'MixedStream', 'HeatUtility', 'PowerUtility', 'Unit', 'units', 'System', 'ureg', 'Quantity', 'Q_', 'units_of_measure']


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
pd.set_option('display.max_columns', 6)
pd.set_option('max_colwidth', 25)

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

from . import exceptions
from . import utils
from . import reports
from . import chemical
from .find import find
from .compound import Compound
from .dissolved_compound import DissolvedCompound
from .chemical import *
from .species import Species
from .stream import Stream
from .mixed_stream import MixedStream
from .heat_utility import HeatUtility
from .power_utility import PowerUtility
from .unit import Unit
from .system import System
from . import units
from .units import *
from .reports import *

__all__.extend(units.__all__)
__all__.extend(reports.__all__)
__all__.extend(chemical.__all__)
SmartBook.Warning = exceptions.DesignWarning
