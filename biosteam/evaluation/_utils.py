# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

import pandas as pd

__all__ = ('var_indices', 'var_columns', 'indices_to_multiindex')


# %% Functions

var_indices = lambda vars: [var.index for var in vars]
var_columns =  lambda vars, names=None: indices_to_multiindex(var_indices(vars), names)
indices_to_multiindex = lambda indices, names=None: pd.MultiIndex.from_tuples(
                                            indices,
                                            names=names or ('Element', 'Variable'),
                                        )
