# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.single_phase_reactor.SinglePhaseReactor

References
----------
.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
    Ng, M. K. Product and Process Design Principles. Wiley 2017.

"""
from .abstract_stirred_tank_reactor import AbstractStirredTankReactor

__all__ = (
    'SinglePhaseReactor', 
)

class SinglePhaseReactor(AbstractStirredTankReactor):
    '''    
    Create a stirred tank reactor which assumes a single outlet stream with 
    no phase equilibrium. 

    The reactor is designed as a pressure vessel with a given aspect ratio and 
    residence time. A pump-heat exchanger recirculation loop can be used to satisfy 
    the duty, if any. By default, a turbine agitator is also included if the 
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
        Power usage of agitator. Defaults to 0.985 [kW / m3] converted from 
        5 hp/1000 gal as in [1]_, for liquid–liquid reaction or extraction.
    vessel_material : 
        Vessel material. Defaults to 'Stainless steel 316'.
    vessel_type : 
        Vessel type. Valid options are 'Horizontal' or 'Vertical'. Defaults to 'Vertical'
    batch :
        Whether to use batch operation mode. If False, operation mode is continuous.
        Defaults to `continuous`.
    tau_0 : 
        Cleaning and unloading time (if batch mode). Defaults to 3 hr.
    N_reactors :
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
    Simulate the liquefaction of corn starch 
    (the enzymatic conversion of starch to glucose):
        
    import biosteam as bst
    from biorefineries import corn
    chemicals = corn.create_chemicals()
    bst.settings.set_thermo(chemicals)
    feed = bst.Stream(
        phase='l', T=360.15, P=1.555e+06, 
        Water=1.115e+05, Ethanol=0.1534, 
        Ash=647.3, Yeast=16.68, CaO=4.73, 
        TriOlein=1754, H2SO4=8.207, NH3=78.84, 
        Starch=2.83e+04, Fiber=5548, 
        SolubleProtein=1899, InsolubleProtein=2318, 
        units='kg/hr'
    )
    R1 = bst.SinglePhaseReactor(
        ins=feed, outs='product', 
        reactions=bst.Reaction('Starch + H2O -> Glucose', 'Starch', 0.95),
        heat_exchanger_configuration='jacketed',
        tau=0.9, V_wf=0.90, V_max=500.,
        kW_per_m3=0.6,
    )
    R1.simulate()
    R1.show('cwt')
    SinglePhaseReactor: R1
    ins...
    [0] feed  
        phase: 'l', T: 360.15 K, P: 1.555e+06 Pa
        composition (%): Water     73.3
                         Ethanol   0.000101
                         Ash       0.426
                         Yeast     0.011
                         CaO       0.00311
                         TriOlein  1.15
                         H2SO4     0.0054
          25.1
                         --------  1.52e+05 kg/hr
    outs...
    [0] product  
        phase: 'l', T: 360.15 K, P: 1.555e+06 Pa
        composition (%): Water     71.4
                         Ethanol   0.000101
                         Glucose   19.6
                         Ash       0.426
                         Yeast     0.011
                         CaO       0.00311
                         TriOlein  1.15
          7.41
                         --------  1.52e+05 kg/hr

    R1.results()
    Single phase reactor                                       Units                   R1
    Electricity         Power                                     kW                 76.9
                        Cost                                  USD/hr                 6.02
    Low pressure steam  Duty                                   kJ/hr             1.83e+06
                        Flow                                 kmol/hr                 47.4
                        Cost                                  USD/hr                 11.3
    Design              Reactor volume                            m3                  142
                        Residence time                            hr                  0.9
                        Vessel type                                              Vertical
                        Length                                    ft                 38.6
                        Diameter                                  ft                 12.9
                        Weight                                    lb             2.38e+05
                        Wall thickness                            in                 1.62
                        Jacketed diameter                                            13.2
                        Vessel material                               Stainless steel 316
    Purchase cost       Vertical pressure vessel (jacketed)      USD             4.09e+05
                        Platform and ladders                     USD              3.6e+04
                        Agitator - Agitator                      USD             5.77e+04
    Total purchase cost                                          USD             5.02e+05
    Utility cost                                              USD/hr                 17.3
    
    '''
    _ins_size_is_fixed = False
    _N_ins = _N_outs = 1

    def _run(self):
        ins = self.ins
        effluent, = self.outs
        if self.P is None: # Constant pressure
            effluent.P = min([i.P for i in ins])
        else:
            effluent.P = self.P
        if self.adiabatic:
            Hnet = sum([i.Hnet for i in ins])
            effluent.mix_from(ins, energy_balance=False)
            self.reactions(effluent)
            effluent.Hnet = Hnet
        else:
            if self.T is None: 
                effluent.mix_from(self.ins) # Constant temperature after mixing
            else:
                effluent.mix_from(self.ins, energy_balance=False)
                effluent.T = self.T
            self.reactions(effluent)