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
    
    {description_doc}
    
    Parameters
    ----------
    {parameters_doc}
        
    Notes
    -----
    {notes_doc}
    
    Examples
    --------
    Simulate the liquefaction of corn starch 
    (the enzymatic conversion of starch to glucose):
        
    >>> import biosteam as bst
    >>> from biorefineries import corn
    >>> chemicals = corn.create_chemicals()
    >>> bst.settings.set_thermo(chemicals)
    >>> feed = bst.Stream(
    ...     phase='l', T=360.15, P=1.555e+06, 
    ...     Water=1.115e+05, Ethanol=0.1534, 
    ...     Ash=647.3, Yeast=16.68, CaO=4.73, 
    ...     TriOlein=1754, H2SO4=8.207, NH3=78.84, 
    ...     Starch=2.83e+04, Fiber=5548, 
    ...     SolubleProtein=1899, InsolubleProtein=2318, 
    ...     units='kg/hr'
    ... )
    >>> R1 = bst.SinglePhaseReactor(
    ...     ins=feed, outs='product', 
    ...     reactions=bst.Reaction('Starch + H2O -> Glucose', 'Starch', 0.95),
    ...     heat_exchanger_configuration='jacketed',
    ...     tau=0.9, V_wf=0.90, V_max=500.,
    ...     kW_per_m3=0.6,
    ... )
    >>> R1.simulate()
    >>> R1.show('cwt')
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
                         ...       25.1
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
                         ...       7.41
                         --------  1.52e+05 kg/hr
    
    >>> R1.results()
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