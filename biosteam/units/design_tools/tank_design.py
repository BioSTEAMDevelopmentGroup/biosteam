# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Design and cost algorithms from ordinary vessels.

References
----------
.. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
    (2009). Techno-economic analysis of a biodiesel production process from
    vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
    https://doi.org/10.1016/j.fuproc.2009.04.017
.. [2] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
    Ng, M. K. Cost Accounting and Capital Cost Estimation.
    In Product and Process Design Principles; Wiley, 2017; pp 426–485.
    
"""
import biosteam as bst
from math import ceil
from thermosteam import settings
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure
from ...utils import ExponentialFunctor

__all__ = ('TankPurchaseCostAlgorithm',
           'field_erected_tank_purchase_cost',
           'compute_number_of_tanks_and_total_purchase_cost',
           'storage_tank_purchase_cost_algorithms',
           'mix_tank_purchase_cost_algorithms')

class TankPurchaseCostAlgorithm:
    r"""
    Create a TankPurchaseCostAlgorithm for vessel costing.
    
    Parameters
    ----------
    f_Cp : function
        Should return the purchase cost given the volume.
    V_min : float
        Minimum volume at which cost is considered accurate.
    V_max : float
        Maximum volume of a vessel.
    V_units : str
        Units of measure for volume.
        
    Attributes
    ----------
    f_Cp : function
        Returns the purchase cost given the volume.
    V_min : float
        Minimum volume at which cost is considered accurate.
    V_max : float
        Maximum volume of a vessel.
    V_units : AbsoluteUnitsOfMeasure
        Units of measure for volume.
    
    Examples
    --------
    Find the number of mixing tanks and the total purchase cost 
    at a volume of 1 m^3 using the purchase cost equation from [1]_:
        
    >>> from biosteam.units.design_tools import TankPurchaseCostAlgorithm
    >>> TankPurchaseCostAlgorithm(lambda V: 12080 * V **0.525,
    ...                             V_min=0.1, V_max=30, V_units='m^3',
    ...                             CE=525.4, material='Stainless steel')
    TankPurchaseCostAlgorithm(f_Cp=<lambda>, V_min=0.1, V_max=30, CE=525.4, material=Stainless steel, V_units=m^3)
    
    
    """
    __slots__ = ('f_Cp', 'V_min', 'V_max',
                 'CE', 'material', 'V_units')

    def __init__(self, f_Cp,  V_min, V_max, V_units, CE, material):
        self.f_Cp = f_Cp
        self.V_min = V_min
        self.V_max = V_max
        self.V_units = AbsoluteUnitsOfMeasure(V_units)
        self.CE = CE
        self.material = material

    def __repr__(self):
        f_name = getattr(self.f_Cp, '__name__', str(self.f_Cp))
        return f"{type(self).__name__}(f_Cp={f_name}, V_min={self.V_min}, V_max={self.V_max}, CE={self.CE}, material={self.material}, V_units={self.V_units})"
    
    
def compute_number_of_tanks_and_total_purchase_cost(total_volume,
                                                    purchase_cost_algorithm):
    """
    Return number of tanks and total purchase cost of all tanks.

    Parameters
    ----------
    total_volume : float
        Total volume required [m^3].
    purchase_cost_algorithm : TankPurchaseCostAlgorithm
        All costing options.
    
    """
    V_total = total_volume
    V_units = purchase_cost_algorithm.V_units
    V_total /= V_units.conversion_factor('m^3')
    V_min = purchase_cost_algorithm.V_min
    if settings.debug and V_total < V_min:
        raise RuntimeError(
             f"volume ({V_total:.5g} {V_units}) is below "
             f"the lower bound ({V_min:.5g} {V_units}) for purchase "
              "cost estimation")
    N = ceil(V_total / purchase_cost_algorithm.V_max)
    if N:
        V = V_total / N
        F_CE = bst.CE / purchase_cost_algorithm.CE
        Cp = N * F_CE * purchase_cost_algorithm.f_Cp(V)
    else:
        Cp = 0.
    return N, Cp
    
def field_erected_tank_purchase_cost(V):
    r"""
    Return the purchase cost [USD] of a single, field-erected vessel assuming
    stainless steel construction material.
    
    Parameters
    ----------
    V : float
        Volume of tank [m^3].
    
    Returns
    -------
    Cp : float
        Purchase cost [USD].
    
    Notes
    -----
    The purchase cost is given by [1]_.
    
    If :math:`V < 2 \cdot 10^3`:
   
        :math:`C_p^{2007} = 32500.0 + 79.35 V`
   
    Otherwise:

        :math:`C_p^{2007} = 125000.0 + 47.1 V`
    
    Examples
    --------
    >>> field_erected_tank_purchase_cost(300)
    112610.0
    
    """
    if V < 2e3:
        Cp = 65000.0 + 158.7 * V
    else:
        Cp = 250000.0 + 94.2 * V
    return Cp

#: Cost algorithms for storage tank vessel types as in [1]_ [2]_.
storage_tank_purchase_cost_algorithms = {
"Field erected": TankPurchaseCostAlgorithm(
    field_erected_tank_purchase_cost,
    V_min=0, V_max=50e3, V_units='m^3',
    CE=525.4, material='Stainless steel'),
"Floating roof": TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=475, n=0.507),
    V_min=3e4, V_max=1e6, V_units='gal',
    CE=567, material='Carbon steel'),
"Cone roof": TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=265, n=0.513),
    V_min=1e4, V_max=1e6, V_units='gal',
    CE=567, material='Carbon steel'),
"Spherical; 0-30 psig": TankPurchaseCostAlgorithm(
    ExponentialFunctor(68, 0.72 ),
    V_min=1e4, V_max=1e6, V_units='gal',
    CE=567, material='Carbon steel'),
"Spherical; 30–200 psig": TankPurchaseCostAlgorithm(
    ExponentialFunctor(53, 0.78),
    V_min=1e4, V_max=7.5e5, V_units='gal',
    CE=567, material='Carbon steel'),
"Gas holder": TankPurchaseCostAlgorithm(
    ExponentialFunctor(3595, 0.43),
    V_min=4e3, V_max=4e5, V_units='ft^3',
    CE=567, material='Carbon steel')
}

#: Cost algorithms for mix tank vessel types as in [2]_.
mix_tank_purchase_cost_algorithms = {
"Conventional": TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=12080, n=0.525),
    V_min=0.1, V_max=30, V_units='m^3',
    CE=525.4, material='Stainless steel')
}