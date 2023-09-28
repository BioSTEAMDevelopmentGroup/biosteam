# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms for the design and purchase cost estimation
of agitators.

References
----------
.. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    Cost Accounting and Capital Cost Estimation (Chapter 16)

"""
from numba import njit
import biosteam as bst
__all__ = (
    'compute_closed_vessel_turbine_purchase_cost',
)

@njit(cache=True)
def _compute_closed_vessel_turbine_purchase_cost(power, CE):
    C_p = 4105 * power ** 0.57
    return CE / 567 * C_p

def compute_closed_vessel_turbine_purchase_cost(power):
    """
    Return the purchase cost [Cp; in USD] of a turbine for closed vessels.
    
    Parameters
    ----------
    power : float
        Power [hp].
    
    Examples
    --------
    >>> compute_closed_vessel_turbine_purchase_cost(power=10)
    15264.97
    
    Notes
    -----
    The purchase cost is given by [1]_. See source code for details.
    The purchase cost is scaled according to BioSTEAM's Chemical
    Plant Cost Index, `biosteam.CE`.
    
    """
    return _compute_closed_vessel_turbine_purchase_cost(power, bst.CE)