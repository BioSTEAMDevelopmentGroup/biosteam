# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 01:40:16 2019

@author: yoelr
"""

from . import biorefinery_tests
from . import binary_distillation_tests

__all__ = (*biorefinery_tests.__all__,
           *binary_distillation_tests.__all__,
           "test_all")

from .biorefinery_tests import *
from .binary_distillation_tests import *

def test_all():
    biorefinery_tests.test_biorefineries()
    binary_distillation_tests.test_all()