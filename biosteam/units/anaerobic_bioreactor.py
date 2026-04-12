# -*- coding: utf-8 -*-
"""
.. contents:: :local:

.. autoclass:: biosteam.units.aerated_bioreactor.AeratedBioreactor
.. autoclass:: biosteam.units.aerated_bioreactor.GasFedBioreactor

References
----------
.. [1] Benz, G. T. Optimize Power Consumption in Aerobic Fermenters. 
    Chem. Eng. Progress 2003, 99 (5), 100–103.

.. [2] Benz, G. T. Bioreactor Design for Chemical Engineers. Chem. Eng.\
    Progress 2011, 21–26.

.. [3] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.

"""
from .abstract_stirred_tank_reactor import (
    AbstractStirredTankReactor, 
    description_doc,
    parameters_doc, 
    notes_doc,
)

__all__ = (
    'AnaerobicBioreactor', 'AnBR',
)

class AnaerobicBioreactor(AbstractStirredTankReactor):
    f"""
    Create an anaerobic bioreactor with a vented stream and an effluent stream .
    
    {description_doc}

    Parameters
    ----------
    {parameters_doc}
        
    Notes
    -----
    {notes_doc}
    
    Examples
    --------
    >>> from biorefineries.cane import create_sugarcane_chemicals
    >>> from biosteam import AnaerobicBioreactor, ReactionSystem, Reaction, Stream, settings
    >>> settings.set_thermo(create_sugarcane_chemicals())
    >>> feed = Stream(Water=1.20e+05,
    ...               Glucose=1.89e+03,
    ...               Sucrose=2.14e+04,
    ...               DryYeast=1.03e+04,
    ...               units='kg/hr',
    ...               T=32+273.15)
    >>> F1 = AnaerobicBioreactor(
    ...     ins=feed, outs=('CO2', 'product'),
    ...     tau=8, heat_exchanger_configuration='jacketed', batch=False,
    ...     reactions=ReactionSystem(
    ...         Reaction('Sucrose + Water -> 2Glucose', 'Sucrose', 1.0),   # Hydrolysis
    ...         Reaction('Glucose -> 2Ethanol + 2CO2',  'Glucose', 0.9),   # Production
    ...         Reaction('Glucose -> Yeast', 'Glucose', 0.70, basis='wt'), # Growth
    ...         basis='mol',
    ...     )
    ... )
    >>> F1.simulate()
    >>> F1.show()
    AnaerobicBioreactor: F1
    ins...
    [0] feed  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.66e+03
                        Glucose  10.5
                        Sucrose  62.5
                        Yeast    456
    outs...
    [0] CO2  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    9.92
                        Ethanol  0.901
                        CO2      244
    [1] product  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): Water    6.59e+03
                        Ethanol  243
                        Glucose  4.07
                        Yeast    532
    
    >>> F1.results()
    Anaerobic bioreactor                                       Units                   F1
    Electricity         Power                                     kW                  301
                        Cost                                  USD/hr                 23.6
    Chilled water       Duty                                   kJ/hr            -1.42e+07
                        Flow                                 kmol/hr              9.5e+03
                        Cost                                  USD/hr                 70.9
    Design              Reactor volume                            m3             1.28e+03
                        Residence time                            hr                    8
                        Vessel type                                              Vertical
                        Length                                    ft                 80.2
                        Diameter                                  ft                 26.7
                        Weight                                    lb             3.35e+05
                        Wall thickness                            in                0.508
                        Jacketed diameter                                            27.1
                        Vessel material                               Stainless steel 316
    Purchase cost       Vertical pressure vessel (jacketed)      USD              5.3e+05
                        Platform and ladders                     USD             1.03e+05
                        Agitator - Agitator                      USD             1.26e+05
    Total purchase cost                                          USD             7.59e+05
    Utility cost                                              USD/hr                 94.5
    
    """
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    T_default = 273.15 + 32 
    P_default = 101325
    jacket_annular_diameter_default = 0.1 # [m]
    V_max_default = 3785 # 1 million gallon is not unheard of for anaerobic bioreactors
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid; reference [1]
    batch_default = True
    
    def _get_duty(self):
        return self.Hnet
    
    @property
    def feed(self):
        return self._ins[0]
    
    @property
    def vent(self):
        return self._outs[0]
        
    def _run_vent(self, vent, effluent):
        vent.phase = 'g'
        vent.receive_vent(effluent, energy_balance=False, ideal=True)
        
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.reactions(effluent)
        for i in self.outs: 
            i.T = self.T
            i.P = self.P
        vent.empty()
        self._run_vent(vent, effluent)

AnBR = AnaerobicBioreactor
