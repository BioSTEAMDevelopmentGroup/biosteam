# -*- coding: utf-8 -*-
"""
.. contents:: :local:

.. autoclass:: biosteam.units.equilibrium_reactor.EquilibriumReactor

"""
import biosteam as bst
from .abstract_stirred_tank_reactor import AbstractStirredTankReactor
from .design_tools.Gibbs_equilibrium_reaction import minimize_Gibbs_free_energy

__all__ = (
    'EquilibriumReactor',
)

class EquilibriumReactor(AbstractStirredTankReactor):
    """
    Create an equilibrium reactor designed as a pressure vessel with a 
    given aspect ratio and residence time. A heat exchanger 
    (e.g., jacket, recirculation loop) can be used to satisfy the duty, if any. 
    By default, a turbine agitator is also included if the 
    power usage, `kW_per_m3`, is positive. A vacuum system is also 
    automatically added if the operating pressure is at a vacuum. 

    Parameters
    ----------
    tau :
        Residence time [hr].
    T : 
        Operating temperature [K].
    P : 
        Operating pressure [Pa].
    V_wf : 
        Fraction of working volume over total volume. Defaults to 0.8.
    length_to_diameter :
        Length to diameter ratio of bioreactor.
    V_max :
        Maximum volume of a reactor [m3]. Defaults to 355.
    kW_per_m3 : 
        Power usage of agitator. Defaults to 0 [kW / m3].
    vessel_material : 
        Vessel material. Defaults to 'Stainless steel 316'.
    vessel_type : 
        Vessel type. Valid options are 'Horizontal' or 'Vertical'. Defaults to 'Vertical'
    batch :
        Whether to use batch operation mode. If False, operation mode is continuous.
        Defaults to `continuous`.
    tau_0 : 
        Cleaning and unloading time (if batch mode). Defaults to 3 hr.
    N :
        Number of reactors.
    heat_exchanger_configuration : 
        What kind of heat exchanger to default to (if any). Valid options include 
        'jacketed', 'recirculation loop', and 'internal coil'. Defaults to 'recirculation loop'.
    dT_hx_loop : 
        Maximum change in temperature for the heat exchanger loop. Defaults to 5 K.
    jacket_annular_diameter :
        Annular diameter of heat exchanger jacket to vessel [m]. Defaults to 0.1 m.
    loading_time :
        Loading time of batch reactor. If not given, it will assume each vessel is constantly
        being filled.
        
    Notes
    -----
    The heat exchanger configuration can be one of the following:

    * 'recirculation loop': 
        The recirculation loop takes into account the required flow rate needed to
        reach the maximum temperature change of the heat exchanger, `dT_hx_loop`. 
        Increasing `dT_hx_loop` decreases the required recirculation flow rate and
        therefore decreases pump costs.
        
        When parallel reactors are required, one recirculation loop (each with a
        pump and heat exchanger) is assumed. Although it is possible to use the
        same recirculation loop for all reactors, this conservative assumption allows
        for each reactor to be operated independently from each other.

    * 'jacketed':
        The jacket does not account for the heat transfer area requirement. 
        It simply assumes that a full jacket can provide the necessary heat transfer
        area to meet the duty requirement. A heuristic annular diameter is assumed
        through `jacket_annular_diameter` (which can be adjusted by the user).
        The temperature at the wall is assumed to be the operating temperature.
        The weight of the jacket is added to the weight of the vessel and the
        cost is compounded together as a jacketed vessel.
        
    * 'internal coil':
        The internal coil is costed as an ordinary helical tube heat exchanger
        with the added assumption that the temperature at the wall is the 
        operating temperature. This method is still not implemented in BioSTEAM
        yet.

    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['N2O4', 'NO2'] pkg='ideal gas')
    >>> feed = bst.Stream(N2O4=100, phase='g', T=273.15 + 25, P=10 * 101325)
    >>> EqR = bst.EquilibriumReactor(
    ...     ins=feed, outs=('gas',), T=feed.T, P=feed.P 
    ... )
    >>> EqR.simulate()
    >>> EqR.show()
    EquilibriumReactor: EqR
    ins...
    [0] feed  
        phase: 'g', T: 298.15 K, P: 1.01325e+06 Pa
        flow: 100 kmol/hr N2O4
    outs...
    [0] gas  
        phase: 'g', T: 298.15 K, P: 1.01325e+06 Pa
        flow (kmol/hr): N2O4  94.1
                        NO2   11.9
    
    >>> EqR.results()
    Equilibrium reactor                                             Units                  EqR
    Low pressure steam       Duty                                   kJ/hr             3.58e+05
                             Flow                                 kmol/hr                 9.24
                             Cost                                  USD/hr                  2.2
    Design                   Reactor volume                            m3                 25.1
                             Number of reactors                                              1
                             Residence time                            hr                  0.1
                             Vessel type                                              Vertical
                             Length                                    ft                 21.6
                             Diameter                                  ft                 7.22
                             Weight                                    lb             3.32e+04
                             Wall thickness                            in                0.705
                             Jacketed diameter                         ft                  7.6
                             Vessel material                               Stainless steel 316
    Purchase cost            Vertical pressure vessel (jacketed)      USD             1.02e+05
                             Platform and ladders                     USD             1.56e+04
    Total purchase cost                                               USD             1.17e+05
    Installed equipment cost                                          USD             4.39e+05
    Utility cost                                                   USD/hr                  2.2
    
    """
    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    T_default = 273.15 + 25
    P_default = 101325
    tau_default = 0.1 # [h]
    jacket_annular_diameter_default = 0.1 # [m]
    heat_exchanger_configuration_default = 'jacketed'
    V_max_default = 500
    kW_per_m3_default = 0
    batch_default = False
    tau_phases_default = 'lg'
    
    def _init(self, phases=None, products=None, **kwargs):
        AbstractStirredTankReactor._init(self, **kwargs)
        self.phases = phases
        self.products = products
    
    def _get_duty(self): return self.Hnet
    
    def _run(self):
        phases = self.phases
        ins = self.ins
        outs = self.outs
        if phases is None:
            phases = []
            for i in ins: phases.extend(i.phases)
            phases = set(phases)
        if len(phases) != len(outs):
            raise RuntimeError(
                'number of phases must be equal to the number of outlets'
            )
        for outlet, phase in zip(outs, phases): outlet.phase = phase
        product = bst.MultiStream.from_streams(outs)
        product.T = self.T
        product.P = self.P
        product.empty()
        outs[0].mix_flows(ins)
        minimize_Gibbs_free_energy(product, self.products)
        
