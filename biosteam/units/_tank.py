# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:47:26 2018

@author: yoelr
"""
from .design_tools import vessel_material_factors
from .._unit import Unit
from ._mixer import Mixer
from math import ceil
import biosteam as bst
from thermosteam import settings
from thermosteam.base import UnitsOfMeasure
from warnings import warn

# %% Cost classes for tanks

ExponentialFunctor = bst.utils.ExponentialFunctor

class VesselPurchaseCostAlgorithm:
    r"""
    Create a VesselPurchaseCostAlgorithm for vessel costing.
    
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
    V_units : UnitsOfMeasure
        Units of measure for volume.
    
    Examples
    --------
    Find the number of mixing tanks and the total purchase cost 
    at a volume of 1 m^3 using the purchase cost equation from [1]_:
        
    >>> from biosteam.units._tank import VesselPurchaseCostAlgorithm
    >>> VesselPurchaseCostAlgorithm(lambda V: 12080 * V **0.525,
    ...                             V_min=0.1, V_max=30, V_units='m^3',
    ...                             CE=525.4, material='Stainless steel')
    VesselPurchaseCostAlgorithm(f_Cp=<lambda>, V_min=0.1, V_max=30, CE=525.4, material=Stainless steel, V_units=m^3)
    
    
    """
    __slots__ = ('f_Cp', 'V_min', 'V_max',
                 'CE', 'material', 'V_units')

    def __init__(self, f_Cp,  V_min, V_max, V_units, CE, material):
        self.f_Cp = f_Cp
        self.V_min = V_min
        self.V_max = V_max
        self.V_units = UnitsOfMeasure(V_units)
        self.CE = CE
        self.material = material

    def __repr__(self):
        return f"{type(self).__name__}(f_Cp={self.f_Cp.__name__}, V_min={self.V_min}, V_max={self.V_max}, CE={self.CE}, material={self.material}, V_units={self.V_units})"


def field_erected_vessel_purchase_cost(V):
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

    Attributes
    ----------
    vessel_type : str
        Vessel type.
    tau : float
        Residence time [hr].
    V_wf : float
        Fraction of working volume over total volume.
    vessel_material : str
        Vessel material.

    Notes
    -----
    The total volume [m^3] is given by:

    :math:`V_{total} = \frac{\tau \cdot Q }{V_{wf}}`

    Where :math:`\tau` is the residence time [hr], :math:`Q` is the flow rate [m^3/hr],
    and :math:`V_{wf}` is the fraction of working volume over total volume.

    The number of tanks is given by:

    :math:`N = Ceil \left( \frac{V_{total}}{V_{max}} \right)`

    Where :math:`V_{max}` is the maximum volume of a tank as given by the cost
    algorithm of the vessel type.

    The volume [m^3] of each tank is given by:

    :math:`V = \frac{V_{total}}{N}`

    The purchase cost will depend on the cost algorithm of the vessel type. 

    Child classes should implement the following class attributes and methods:
    
    purchase_cost_algorithms : dict[str: VesselPurchaseCostAlgorithm]
        All purchase cost algorithms available for vessel types.

    """
    _units = {'Total volume': 'm^3'}
    _N_outs = 1

    def __init_subclass__(cls, isabstract=False):
        if not isabstract:
            if not hasattr(cls, 'purchase_cost_algorithms'):
                raise NotImplementedError("Tank subclass must implement "
                                          "a 'purchase_cost_algorithms' dictionary")
        super().__init_subclass__(isabstract)

    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, material):
        self._vessel_material = material
        default_material = self.purchase_cost_algorithm.material
        if material == default_material:
            self._F_M = 1.0
        else:
            try:
                F_M_new = vessel_material_factors[material]
            except:
                raise ValueError("no material factor available for "
                                f"vessel construction material '{material}';"
                                 "only the following materials are "
                                f"available: {', '.join(vessel_material_factors)}")
            try:
                F_M_old = vessel_material_factors[default_material]
            except KeyError:
                raise ValueError(f"only '{default_material}' is a valid "
                                    "material for given vessel type")
            self._F_M = F_M_new / F_M_old        
    
    @property
    def vessel_type(self):
        return self._vessel_type
    @vessel_type.setter
    def vessel_type(self, vessel_type):
        if vessel_type in self.purchase_cost_algorithms:
            self._vessel_type = vessel_type
            self.purchase_cost_algorithm = self.purchase_cost_algorithms[vessel_type]
        else:
            raise ValueError(f"vessel type '{vessel_type}'"
                             "is not avaiable; only the following vessel "
                             "types are available: "
                            f"{', '.join(self.purchase_cost_algorithms)}")
    
    def _design(self):
        design_results = self.design_results
        design_results['Residence time'] = tau = self.tau
        design_results['Total volume'] = tau * self.F_vol_out / self.V_wf

    def _cost(self):
        V_total = self.design_results['Total volume']
        Cp_algorithm = self.purchase_cost_algorithm
        V_units = Cp_algorithm.V_units
        V_total /= V_units.conversion_factor('m^3')
        if settings.debug and V_total < Cp_algorithm.V_min:
            warn(f"volume of {self} ({V_total:.5g} {V_units}) is below "
                 f"the lower bound ({self.V_min:.5g} {V_units}) for purchase "
                  "cost estimation", RuntimeWarning)
        N = ceil(V_total / Cp_algorithm.V_max)
        V = V_total / N
        Cp = Cp_algorithm.f_Cp(V)
        F_CE = bst.CE / Cp_algorithm.CE
        self.design_results['Number of tanks'] = N
        self.purchase_costs['Tanks'] = N * self._F_M * F_CE * Cp


# %% Storage tank purchase costs    


class StorageTank(Tank):
    r"""
    Create a storage tank with a purchase cost based on volume as given by residence time.

    Parameters
    ----------
    ins : stream
        Inlet.
    outs : stream
        Outlet.
    vessel_type : str, optional
        Vessel type. Defaults to 'Field erected'.
    tau=672 : float
        Residence time [hr].
    V_wf=1 : float
        Fraction of working volume over total volume.
    vessel_material : str, optional
        Vessel material. Defaults to 'Stainless steel'.

    Notes
    -----
    For a detailed discussion on the design and cost algorithm,
    please read the :doc:`Tank` documentation.
    
    References to the purchase cost algorithms for each available
    vessel type are as follows:
        
    * Field erected: [1]_
    * Floating roof: [2]_
    * Cone roof: [2]_
    
    Examples
    --------
    Create a carbon steel, floating roof storage tank for the storage of bioethanol:
    
    >>> import thermosteam as tmo
    >>> from biosteam import units
    >>> tmo.settings.set_thermo(tmo.Chemicals(['Ethanol']))
    >>> feed = tmo.Stream('feed', Ethanol=23e3, units='kg/hr')
    >>> effluent = tmo.Stream('effluent')
    >>> T1 = units.StorageTank('T1', ins=feed, outs=effluent,
    ...                        tau=7*24, # In hours
    ...                        vessel_type='Floating roof',
    ...                        vessel_material='Carbon steel')
    >>> T1.simulate()
    >>> T1.show(flow='kg/hr')
    StorageTank: T1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Ethanol  2.3e+04
    outs...
    [0] effluent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Ethanol  2.3e+04
    
    >>> T1.results()
    Tank                                  Units       T1
    Design              Residence time               168
                        Total volume        m^3 4.92e+03
                        Number of tanks                2
    Purchase cost       Tanks               USD 8.42e+05
    Total purchase cost                     USD 8.42e+05
    Utility cost                         USD/hr        0
    
    
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
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type='Field erected', tau=4*7*24,
                  V_wf=1.0, vessel_material='Stainless steel'):
        Unit.__init__(self, ID, ins, outs, thermo)

        # [str] Vessel type.
        self.vessel_type = vessel_type

        #: [float] Residence time in hours.
        self.tau = tau

        #: [float] Fraction of working volume to total volume.
        self.V_wf = V_wf

        #: [str] Vessel construction material.
        self.vessel_material = vessel_material
    
    #: dict[str: VesselPurchaseCostAlgorithm] All cost algorithms available for vessel types.
    purchase_cost_algorithms = {
    "Field erected": VesselPurchaseCostAlgorithm(
        field_erected_vessel_purchase_cost,
        V_min=0, V_max=50e3, V_units='m^3',
        CE=525.4, material='Stainless steel'),
    "Floating roof": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(A=475, n=0.507),
        V_min=3e4, V_max=1e6, V_units='gal',
        CE=567, material='Carbon steel'),
    "Cone roof": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(A=265, n=0.513),
        V_min=1e4, V_max=1e6, V_units='gal',
        CE=567, material='Carbon steel'),
    "Spherical; 0-30 psig": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(68, 0.72 ),
        V_min=1e4, V_max=1e6, V_units='gal',
        CE=567, material='Carbon steel'),
    "Spherical; 30–200 psig": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(53, 0.78),
        V_min=1e4, V_max=7.5e5, V_units='gal',
        CE=567, material='Carbon steel'),
    "Gas holder": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(3595, 0.43),
        V_min=4e3, V_max=4e5, V_units='ft^3',
        CE=567, material='Carbon steel')
    }


class MixTank(Tank):
    """
    Create a mixing tank with volume based on residence time.

    Parameters
    ----------
    ins : streams
        Inlet fluids to be mixed.
    outs : stream
        Outlet.
    vessel_type : str
        Vessel type.
    tau=1 : float
        Residence time [hr].
    V_wf=0.8 : float
        Fraction of working volume over total volume.
    vessel_material : str, optional
        Vessel material. Defaults to 'Stainless steel'.
    kW_per_m3=0.0985 : float
        Electricity requirement per unit volume [kW/m^3].
    
    Notes
    -----
    For a detailed discussion on the design and cost algorithm,
    please read the :doc:`Tank` documentation.
    
    The purchase cost algorithm is based on [1]_.

    References
    ----------
    .. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
        (2009). Techno-economic analysis of a biodiesel production process from
        vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
        https://doi.org/10.1016/j.fuproc.2009.04.017

    """
    _N_ins = 2
    _run = Mixer._run
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type="Conventional", tau=1,
                  V_wf=0.8, vessel_material='Stainless steel',
                  kW_per_m3=0.0985):
        Unit.__init__(self, ID, ins, outs, thermo)
        
        # [float] Electricity requirement per unit volume [kW/m^3].
        self.kW_per_m3 = kW_per_m3
        
        # [str] Vessel type.
        self.vessel_type = vessel_type

        #: [float] Residence time in hours.
        self.tau = tau

        #: [float] Fraction of working volume to total volume.
        self.V_wf = V_wf

        #: [str] Vessel construction material.
        self.vessel_material = vessel_material
    
    #: dict[str: VesselPurchaseCostAlgorithm] All cost algorithms available for vessel types.
    purchase_cost_algorithms = {
    "Conventional": VesselPurchaseCostAlgorithm(
        ExponentialFunctor(A=12080, n=0.525),
        V_min=0.1, V_max=30, V_units='m^3',
        CE=525.4, material='Stainless steel')
    }

    def _cost(self):
        super()._cost()
        self.power_utility(self.kW_per_m3 * self.design_results['Total volume'])


del ExponentialFunctor