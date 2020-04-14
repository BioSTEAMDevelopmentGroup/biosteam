# -*- coding: utf-8 -*-
from thermosteam import exceptions
from thermosteam.exceptions import *

__all__ = ('DesignError', *exceptions.__all__)
del exceptions

# %% Biosteam errors

class DesignError(RuntimeError):
    """RuntimeError regarding unit design."""

