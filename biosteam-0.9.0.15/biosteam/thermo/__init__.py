# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

from . import acentric
from . import activity
from . import chemical
from . import combustion
from . import critical
from . import coolprop
from . import dipole
from . import dippr
from . import datasheet
from . import electrochem
from . import elements
from . import environment
from . import eos
from . import eos_mix
from . import heat_capacity
from . import identifiers
from . import joback
from . import law
from . import lennard_jones
from . import miscdata
from . import mixture
from . import permittivity
from . import phase_change
from . import property_package
from . import reaction
from . import refractivity
from . import safety
from . import solubility
from . import stream
from . import interface
from . import thermal_conductivity
from . import triple
from . import unifac
from . import utils
from . import vapor_pressure
from . import virial
from . import viscosity
from . import volume

__all__ = ['activity', 'chemical', 'combustion', 'critical',
           'dipole', 'electrochem', 'elements', 'environment', 'eos', 'eos_mix',
           'heat_capacity',  'identifiers', 'joback', 'law', 'lennard_jones',
           'miscdata', 'permittivity', 'phase_change', 'property_package', 'reaction',
           'refractivity', 'safety', 'solubility', 'interface',
           'thermal_conductivity', 'triple', 'utils', 'vapor_pressure', 'virial', 'viscosity', 
           'volume', 'acentric', 'coolprop', 'datasheet', 'dippr', 'unifac', 'stream', 'mixture']

from lazypkg import LazyPkg
LazyPkg(__name__, ())

__version__ = '0.1.39'
