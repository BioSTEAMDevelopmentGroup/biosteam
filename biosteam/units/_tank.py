# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:47:26 2018

@author: yoelr
"""
from .designtools import vessel_material_factors
from .._unit import Unit
from .._exceptions import DesignError
from ._mixer import Mixer
from .decorators import cost
from math import ceil
import numpy as np
import biosteam as bst
from thermosteam.base import UnitsOfMeasure
from warnings import warn

# %% Cost classes for tanks

class ExponentialFunction:
    __slots__ = ('base', 'n')

    def __init__(self, base, n):
        self.base = base
        self.n = n

    def __call__(self, S):
        return self.base * S ** self.n

    def __repr__(self):
        return f"{type(self).__name__}(base={self.base}, n={self.n})"


class VesselCostResults:
    __slots__ = ('N', 'Cp', 'kW')
    def __init__(self, N, Cp, kW):
        self.N = N
        self.Cp = Cp
        self.kW = kW

    def __repr__(self):
        return f"{type(self).__name__}(N={self.N}, Cp={self.Cp}, kW={self.kW})"


class VesselCostAlgorithm:
    """
    Create a VesselCostAlgorithm for vessel costing.
    
    Parameters
    ----------
    f_Cp : function
        Should return the purchase cost given the volume.
    kW : float
        Electricity rate [kW] per unit volume.
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
    kW : float
        Electricity rate [kW] per unit volume.
    V_min : float
        Minimum volume at which cost is considered accurate.
    V_max : float
        Maximum volume of a vessel.
    V_units : UnitsOfMeasure
        Units of measure for volume.
    
    Examples
    --------
    Find the number of mixing tanks, the total purchase cost, and the total
    electricity rate given a volume of 1 m^3 according to [1]_:
        
    >>> from biosteam.units._tank import VesselCostAlgorithm
    >>> cost_algorithm = VesselCostAlgorithm(lambda V: 12080 * V **0.525,
                                             V_min=0.1, V_max=30, V_units='m^3',
                                             CE=525.4, kW=0.0985, 
                                             material='Stainless steel')
    >>> results = cost_algorithm(V_total=1., units='m^3', material='Stainless steel')
    VesselCostResults(N=1, Cp=13047.963456414163, kW=0.0985)
    
    """
    __slots__ = ('f_Cp', 'kW', 'V_min', 'V_max',
                 'CE', 'material', 'V_units')
    vessel_material_factors = vessel_material_factors

    def __init__(self, f_Cp, kW, V_min, V_max, V_units, CE, material):
        self.f_Cp = f_Cp
        self.kW = kW
        self.V_min = V_min
        self.V_max = V_max
        self.V_units = UnitsOfMeasure(V_units)
        self.CE = CE
        self.material = material
        

    def __call__(self, V_total, units, material, source="tank"):
        V_units = self.V_units
        V_total /= V_units.conversion_factor(units)
        if V_total < self.V_min:
            warn(f"volume of {repr(source)} ({V_total:.5g} {V_units}) is below "
                 f"the lower bound ({self.V_min:.5g} {V_units})",
                 RuntimeWarning)
        kW = self.kW * V_total
        N = ceil(V_total / self.V_max)
        V = V_total / N
        Cp = self.f_Cp(V)
        material_factors = self.vessel_material_factors
        default_material = self.material
        if default_material == material:
            F_M = 1.0
        else:
            F_M = material_factors[material] / material_factors[default_material]
        F_CE = bst.CE / self.CE
        return VesselCostResults(N, N * F_M * F_CE * Cp, kW)

    def __repr__(self):
        return f"{type(self).__name__}(f_Cp={self.f_Cp}, kW={self.kW}, V_min={self.V_min}, V_max={self.V_max}, CE={self.CE}, material={self.material}, V_units={self.V_units})"


def field_erected_vessel_purchase_cost(V):
    """
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

    .. math::

        If V < 2 \cdot 10^3:
   
            C_p^{2007} &= 32500.0 + 79.35 V
   
        Otherwise:
    
            C_p^{2007} &= 125000.0 + 47.1 V

    Examples
    --------
    >>> from biosteam.units._tank import field_erected_vessel_purchase_cost
    >>> field_erected_vessel_purchase_cost(300)
    112610.0

    References
    ----------
    .. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
        (2009). Techno-economic analysis of a biodiesel production process from
        vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
        https://doi.org/10.1016/j.fuproc.2009.04.017

    """
    if V < 2e3:
        Cp = 65000.0 + 158.7 * V
    else:
        Cp = 250000.0 + 94.2 * V
    return Cp



# %%

class Tank(Unit, isabstract=True):
    r"""
    Abstract class for tanks.

    Parameters
    ----------
    tau : float
        Residence time [hr].
    V_wf : float
        Fraction of working volume over total volume.
    material : str
        Vessel material.

    Notes
    -----
    The total volume [m^3] is given by:

    :math:`V_{total} = \frac{\tau \cdot Q }{V_f}`

    Where :math:`\tau` is the residence time [hr], :math:`Q` is the flow rate [m^3/hr],
    and :math:`V_f` is the fraction of working volume over total volume.

    The number of tanks is given by:

    :math:`N = Ceil \left( \frac{V_{total}}{V_{max}} \right)`

    Where :math:`V_{max}` is the maximum volume of a tank as given by the cost
    algorithm of the vessel type.

    The volume [m^3] of each tank is given by:

    :math:`V = \frac{V_{total}}{N}`

    The purchase cost will depend on the cost algorithm of the vessel type. 

    Child classes should implement the following class attributes and methods:

    defaults : dict[str: value]
        Default values for `tau`, `V_wf`, and `material`.
    
    cost_algorithms : dict[str: VesselCostAlgorithm]
        All cost algorithms available for vessel types.

    Attributes
    ----------
    V_wf : float
        Fraction of working volume over total volume [m^3].
    tau : float
        Residence time [hr].
    material : str
        Vessel material.

    """
    _units = {'Total volume': 'm^3'}
    _N_outs = 1

    def __init_subclass__(cls, isabstract=False):
        if not isabstract:
            try:
                defaults = cls.defaults
            except:
                raise NotImplementedError("Tank subclass must implement a "
                                          "'defaults' dictionary")
            for key in ('vessel_type', 'tau', 'V_wf', 'material'):
                if key not in defaults:
                    raise NotImplementedError("'defaults' dictionary must "
                                             f"include a '{key}' item")
            if not hasattr(cls, 'cost_algorithms'):
                raise NotImplementedError("Tank subclass must implement "
                                          "a 'cost_algorithms' dictionary")
        super().__init_subclass__(isabstract)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type=None, tau=None, V_wf=None, material=None):
        Unit.__init__(self, ID, ins, outs, thermo)

        defaults = self.defaults
        self.vessel_type = vessel_type or defaults['vessel_type']

        #: [float] Residence time in hours.
        self.tau = tau or defaults['tau']

        #: [float] Fraction of working volume to total volume.
        self.V_wf = V_wf or defaults['V_wf']

        #: [str] Vessel construction material.
        self.material = material or defaults['material']

    @property
    def cost_algorithm(self):
        """[VesselCostAlgorithm] The algorithm and parameters being used to
        calculate the purchase cost."""
        return self.cost_algorithms[self.vessel_type]

    def _design(self):
        design_results = self.design_results
        design_results['Residence time'] = tau = self.tau
        design_results['Total volume'] = tau * self.F_vol_out / self.V_wf

    def _cost(self):
        V_total = self.design_results['Total volume']
        results = self.cost_algorithm(V_total, 'm^3', self.material, self)
        self.design_results['Number of tanks'] = results.N
        self.purchase_costs['Tanks'] = results.Cp
        self.power_utility(results.kW)


# %% Storage tank purchase costs    


class StorageTank(Tank):
    r"""
    Create a storage tank with a purchase cost based on volume as given by residence time.

    Parameters
    ----------
    ins : 
        [0] Feed stream.
    outs : 
        [0] Effluent stream.

    Notes
    -----
    For a detailed discussion on the design and cost algorithm,
    please read the :doc:`Tank` documentation.
    
    References for the purchase cost algorithms at a given vessel type are as follows:
        
    * Field erected: [1]_
    * Floating roof: [2]_
    * Cone roof: [2]_
    
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
    #: dict[str: VesselCostAlgorithm] All cost algorithms available for vessel types.
    cost_algorithms = {
    "Field erected": VesselCostAlgorithm(
        field_erected_vessel_purchase_cost, kW=0,
        V_min=0, V_max=50e3, V_units='m^3', CE=525.4,
        material='Stainless steel'),
    "Floating roof": VesselCostAlgorithm(
        ExponentialFunction(base=475, n=0.507), kW=0,
        V_min=3e4, V_max=1e5, V_units='gal', CE=567,
        material='Carbon steel'),
    "Cone roof": VesselCostAlgorithm(
        ExponentialFunction(base=265, n=0.513), kW=0,
        V_min=1e4, V_max=1e5, V_units='gal', CE=567,
        material='Carbon steel'),
    }

    #: [dict] All default values.
    defaults = dict(
        vessel_type="Field erected", 
        tau=4*7*24, 
        V_wf=1.0, 
        material = 'Stainless steel'
    )



class MixTank(Tank):
    """
    Create a mixing tank with volume based on residence time.

    Parameters
    ----------
    ins : 
        [:] Feed streams.
    outs : 
        [0] Effluent stream.

    Notes
    -----
    For a detailed discussion on the design and cost algorithm,
    please read the :doc:`Tank` documentation.
    
    The purchase cost algorithms are based on [1]_.

    References
    ----------
    .. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
        (2009). Techno-economic analysis of a biodiesel production process from
        vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
        https://doi.org/10.1016/j.fuproc.2009.04.017

    """
    _N_ins = 2
    _run = Mixer._run

    #: dict[str: VesselCostAlgorithm] All cost algorithms available for vessel types.
    cost_algorithms = {
    "Conventional": VesselCostAlgorithm(
        ExponentialFunction(base=12080, n=0.525),
        V_min=0.1, V_max=30, V_units='m^3',
        CE=525.4, kW=0.0985, 
        material='Stainless steel')
    }

    #: [dict] All default values.
    defaults = dict(
        vessel_type="Conventional", 
        tau=1, 
        V_wf=0.8, 
        material='Stainless steel'
    )


