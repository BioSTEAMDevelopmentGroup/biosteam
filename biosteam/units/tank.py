# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for tanks.

.. contents:: :local:
    
Unit operations
---------------
.. autoclass:: biosteam.units.tank.Tank
.. autoclass:: biosteam.units.tank.MixTank
.. autoclass:: biosteam.units.tank.StorageTank

Factories
---------
.. autofunction:: biosteam.units.tank.tank_factory

"""
from .design_tools.specification_factors import vessel_material_factors
from .design_tools.tank_design import (
    TankPurchaseCostAlgorithm,
    compute_number_of_tanks_and_total_purchase_cost,
    storage_tank_purchase_cost_algorithms,
    mix_tank_purchase_cost_algorithms)
from ..utils import ExponentialFunctor
from .._unit import Unit
from .mixing import Mixer

__all__ = ('Tank', 'MixTank', 'StorageTank', 'tank_factory')

# %% Factory functions

def tank_factory(name, *, CE, cost, S, tau, n=0.6, kW_per_m3=0., V_wf=0.9, 
                 V_min=0., V_max=1e6, V_units='m3', material='Carbon steel', 
                 vessel_type='Tank', BM=None, module=None, mixing=False):
    """
    Return a Tank sublcass.

    Parameters
    ----------
    name : str
        Name of subclass.
    CE : float
        Chemical plant cost index.
    cost : float
        Cost of tank.
    S : float
        Size of tank.
    tau : float
        Residence time [hr].
    n : float, optional
        Scale up factor. The default is 0.6.
    kW_per_m3 : float, optional
        Electricity consumption per volume [kW/m3]. The default is 0..
    V_wf : float, optional
        Fraction of working volume. The default is 0.9.
    V_min : float, optional
        Minimum tank volume. The default is 0..
    V_max : float, optional
        Maximum tank volume. The default is 1e6.
    V_units : str, optional
        Volume units of measure. The default is 'm3'.
    material : str, optional
        Vessel material. The default is 'Carbon steel'.
    vessel_type : str, optional
        Name of vessel type. The default is 'Tank'.
    BM : float, optional
        Bare module factor. The default is 2.3 for all tanks.
    module : str, optional
        Module to implement class.
    mixing: bool, optional
        Whether multiple streams are mixed in the tank.
        
    Examples
    --------
    >>> import biosteam as bst
    >>> Corn = bst.Chemical('Corn', search_db=False, default=True, phase='s')
    >>> bst.settings.set_thermo([Corn])
    >>> CornStorage = bst.tank_factory('CornStorage', 
    ...     CE=525, cost=979300., S=185400, tau=259.2, n=1.0, V_wf=0.9, V_max=3e5, 
    ...     V_units='m3'
    ... )
    >>> corn = bst.Stream('corn', Corn=46350.72444)
    >>> T101 = CornStorage('T101', corn, 'outlet')
    >>> T101.simulate()
    >>> T101.show()
    CornStorage: T101
    ins...
    [0] corn
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Corn  4.64e+04
    outs...
    [0] outlet
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Corn  4.64e+04
    >>> T101.results()
    Corn storage                          Units     T101
    Design              Residence time       hr      259
                        Total volume        m^3 1.33e+04
                        Number of tanks                1
    Purchase cost       Tanks               USD 7.62e+04
    Total purchase cost                     USD 7.62e+04
    Utility cost                         USD/hr        0

    """
    dct = {
        'purchase_cost_algorithms': {
            vessel_type: TankPurchaseCostAlgorithm(
                ExponentialFunctor(cost/S**n, n), V_min, V_max, V_units, CE, material
            )
        },
        '_default_vessel_type': vessel_type,
        '_default_tau': tau,
        '_default_V_wf': V_wf,
        '_default_vessel_material': material,
        '_default_kW_per_m3': kW_per_m3,
    }
    if BM: 
        dct['_F_BM_default'] = {'Tanks': BM}
    if mixing:
        dct['_ins_size_is_fixed'] = False
        dct['_N_ins'] = 2
        dct['_run'] = Mixer._run
    cls = type(name, (Tank,), dct)
    if module: cls.__module__ = module
    return cls

# %% Abstract tank class

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
    kW_per_m3 : float
        Electricity requirement per unit volume [kW/m^3].

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
    _default_kW_per_m3 = 0.
    _units = {'Total volume': 'm^3',
              'Residence time': 'hr'}
    _F_BM_default = {'Tanks': 2.3}
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type=None, tau=None, V_wf=None, 
                  vessel_material=None, kW_per_m3=0.):
        Unit.__init__(self, ID, ins, outs, thermo)

        # [str] Vessel type.
        self.vessel_type = vessel_type or self._default_vessel_type

        #: [float] Residence time in hours.
        self.tau = tau or self._default_tau

        #: [float] Fraction of working volume to total volume.
        self.V_wf = V_wf or self._default_V_wf

        #: [str] Vessel construction material.
        self.vessel_material = vessel_material or self._default_vessel_material
        
        # [float] Electricity requirement per unit volume [kW/m^3].
        self.kW_per_m3 = kW_per_m3 or self._default_kW_per_m3

    def __init_subclass__(cls, isabstract=False):
        if not isabstract:
            hasfield = hasattr
            if not hasfield(cls, 'purchase_cost_algorithms'):
                raise NotImplementedError("Tank subclass must implement "
                                          "a 'purchase_cost_algorithms' dictionary")
            attributes = ('_default_vessel_type', '_default_tau', 
                          '_default_V_wf', '_default_vessel_material')
            for i in attributes:
                if not hasfield(cls, i):
                    raise NotImplementedError("Tank subclass must implement "
                                              "a '{i}' attribute")
        super().__init_subclass__(isabstract)

    @property
    def vessel_material(self):
        return self._vessel_material
    @vessel_material.setter
    def vessel_material(self, material):
        self._vessel_material = material
        default_material = self.purchase_cost_algorithm.material
        if material == default_material:
            self.F_M['Tanks'] = vessel_material_factors.get(default_material, 1.)
        else:
            try:
                self.F_M['Tanks'] = vessel_material_factors[material]
            except:
                raise ValueError("no material factor available for "
                                f"vessel construction material '{material}';"
                                 "only the following materials are "
                                f"available: {', '.join(vessel_material_factors)}")
    
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
        design_results = self.design_results
        V = design_results['Total volume']
        design_results['Number of tanks'], Cp = compute_number_of_tanks_and_total_purchase_cost(
            V, self.purchase_cost_algorithm
        )
        self.baseline_purchase_costs['Tanks']  = Cp / vessel_material_factors.get(self._vessel_material, 1.)
        self.power_utility.consumption = self.kW_per_m3 * V


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
    tau : float, optional
        Residence time [hr]. Defaults to 672.
    V_wf : float, optional
        Fraction of working volume over total volume. Defaults to 1.
    vessel_material : str, optional
        Vessel material. Defaults to 'Stainless steel'.
    kW_per_m3 : float, optional
        Electricity requirement per unit volume [kW/m^3]. Defaults to 0.

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
    >>> settings.set_thermo(['Ethanol'], cache=True)
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
    Design              Residence time       hr      168
                        Total volume        m^3 4.92e+03
                        Number of tanks                2
    Purchase cost       Tanks               USD 8.41e+05
    Total purchase cost                     USD 8.41e+05
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
    _outs_size_is_fixed = _ins_size_is_fixed = True
    _default_vessel_type = 'Field erected'
    _default_tau = 4*7*24
    _default_V_wf = 1.0
    _default_vessel_material = 'Stainless steel'
    _default_kW_per_m3 = 0.
    
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
    vessel_type : str, optional
        Vessel type. Defaults to 'Conventional'.
    tau : float, optional
        Residence time [hr]. Defaults to 1.
    V_wf : float, optional
        Fraction of working volume over total volume. Defaults to 0.8.
    vessel_material : str, optional
        Vessel material. Defaults to 'Stainless steel'.
    kW_per_m3 : float, optional
        Electricity requirement per unit volume [kW/m^3]. Defaults to 0.0985.
    
    Notes
    -----
    For a detailed discussion on the design and cost algorithm,
    please read the :doc:`Tank` documentation.
    
    The purchase cost algorithm is based on [1]_.
    The electricity rate is based on [2]_.

    References
    ----------
    .. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
        (2009). Techno-economic analysis of a biodiesel production process from
        vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
        https://doi.org/10.1016/j.fuproc.2009.04.017

    .. [2] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)

    """
    _N_ins = 2
    _ins_size_is_fixed = False
    _run = Mixer._run
    _default_vessel_type = 'Conventional'
    _default_tau = 1
    _default_V_wf = 0.8
    _default_vessel_material = 'Stainless steel'
    _default_kW_per_m3 = 0.0985
    
    #: dict[str: VesselPurchaseCostAlgorithm] All cost algorithms available for vessel types.
    purchase_cost_algorithms = mix_tank_purchase_cost_algorithms


MixTank._graphics.edge_in *= 3
MixTank._graphics.edge_out *= 3
