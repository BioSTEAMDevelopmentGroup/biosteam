# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 01:40:16 2019

@author: yoelr
"""

from . import biorefinery_tests
from . import unit_operation_tests
from . import utility_tests

__all__ = (*biorefinery_tests.__all__,
           "test_biosteam")

from .biorefinery_tests import *
from .unit_operation_tests import *
from .utility_tests import *

def test_biosteam():
    utility_tests.test_utilities()
    unit_operation_tests.test_unit_operations()
    biorefinery_tests.test_biorefineries()