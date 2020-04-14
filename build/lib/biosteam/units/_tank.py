# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:47:26 2018

@author: yoelr
"""
from .design_tools.specification_factors import vessel_material_factors
from .design_tools.tank_design import (
    compute_number_of_tanks_and_total_purchase_cost,
    storage_tank_purchase_cost_algorithms,
    mix_tank_purchase_cost_algorithms)
from ..utils import ExponentialFunctor
from .._unit import Unit
from ._mixer import Mixer
import biosteam as bst

from warnings import warn

__all__ = ('Tank', 'MixTank', 'StorageTank')

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
    The total volume [:math:`m^3`] is given by:

    :math:`V_{total} = \frac{\tau \cdot Q }{V_{wf}}`

    Where :math:`\tau` is the residence time [:math:`hr`], 
    :math:`Q` is the flow rate [:math:`\frac{m^3}{hr}`],
    and :math:`V_{wf}` is the fraction of working volume over total volume.

    The number of tanks is given by:

    :math:`N = Ceil \left( \frac{V_{total}}{V_{max}} \right)`

    Where :math:`V_{max}` is the maximum volume of a tank depends on the 
    vessel type.

    The volume [:math:`m^3`] of each tank is given by:

    :math:`V = \frac{V_{total}}{N}`

    The purchase cost will depend on the vessel type. 

    Child classes should implement the following class attributes and methods:
    
    purchase_cost_algorithms : dict[str: VesselPurchaseCostAlgorithm]
        All purchase cost algorithms options for selected vessel types.

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
        N, Cp = compute_number_of_tanks_and_total_purchase_cost(
            self.design_results['Total volume'],
            self.purchase_cost_algorithm,
            self._F_M)
        self.design_results['Number of tanks'] = N
        self.purchase_costs['Tanks'] = Cp


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
    
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(['Ethanol'])
    >>> feed = Stream('feed', Ethanol=23e3, units='kg/hr')
    >>> effluent = Stream('effluent')
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
    Storage tank                          Units       T1
    Design              Residence time               168
                        Total volume        m^3 5.19e+03
                        Number of tanks                2
    Purchase cost       Tanks               USD 8.65e+05
    Total purchase cost                     USD 8.65e+05
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
    purchase_cost_algorithms = storage_tank_purchase_cost_algorithms


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
    purchase_cost_algorithms = mix_tank_purchase_cost_algorithms

    def _cost(self):
        super()._cost()
        self.power_utility(self.kW_per_m3 * self.design_results['Total volume'])

MixTank._graphics.edge_in *= 3
MixTank._graphics.edge_out *= 3
