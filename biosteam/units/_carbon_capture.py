# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com> and Yalin Li <zoe.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from .._unit import Unit
from .decorators import cost
from .design_tools.cost_index import CEPCI_by_year

__all__ = ('AmineAbsorption', 'CO2Compression')

@cost(basis='Total flow', ID='Absorber', units='kmol/hr',
      cost=4.81e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.3)
@cost(basis='Total flow', ID='Stripper', units='kmol/hr',
      cost=4e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.3)
@cost(basis='CO2 flow', ID='Pumps', units='kmol/hr',
      # 55595.96 for 613 metric tonne/hr CO2
      kW=55595.96/(613*(1000/44))*(24123*0.1186),
      cost=0.42e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.3)
@cost(basis='Total flow', ID='Condenser', units='kmol/hr',
      cost=0.27e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=4.17)
@cost(basis='Total flow', ID='Reboiler', units='kmol/hr',
      cost=0.53e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Cross heat exchanger', units='kmol/hr',
      cost=2.28e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Cooler', units='kmol/hr',
      cost=0.09e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.17)
@cost(basis='Total flow', ID='Makeup tank', units='kmol/hr',
      cost=0.23e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=2.3)
class AmineAbsorption(Unit):
    '''
    Create an AmineAbsorption unit for capture of CO2 in the flue gas using
    30 wt% aqueous monoethanolamine (MEA). Capital cost and duty basis are
    from the conventional configuration as detailed in [1]_.
    Cost is extrapolated via the 6/10th rule [2]_. Pump power usage is
    based on [2]_. Bare module factors are based on similar units in BioSTEAM.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] Flue gas containing CO2
        * [1] Makeup MEA (neat), updated by the unit
        * [2] Makeup water, updated by the unit
    outs : stream sequence
        * [0] CO2-stripped vent
        * [1] Concentrated CO2
    CO2_recovery :
        Percentage of CO2 that can be captured.
    MEA_to_CO2 :
        Net usage of MEA (kg pure MEA/metric tonne CO2 captured).
        The default is 1.5 based on [1]_ and [3]_.
    heat_ratio :
        Unit duty in kJ/kg CO2.
    
    Examples
    --------
    
    >>> import biosteam as bst
    >>> import thermosteam as tmo
    >>> MEA = tmo.Chemical('MEA', search_ID='141-43-5', phase='l')
    >>> chems = tmo.Chemicals(('CO2', 'O2', 'Water', 'N2', MEA))
    >>> tmo.settings.set_thermo(chems)
    >>> flue_gas = tmo.Stream('flue_gas',
    ...                       CO2=2895,
    ...                       N2=17609,
    ...                       O2=3618,
    ...                       units='kmol/hr',
    ...                       T=48+273.15)
    >>> U1 = bst.units.AmineAbsorption('U1',
    ...                                ins=(flue_gas, 'makeup_MEA', 'makeup_water'),
    ...                                outs=('vent', 'CO2'))
    >>> U1.simulate()
    >>> U1.show()
    AmineAbsorption: U1
    ins...
    [0] flue_gas
        phase: 'l', T: 321.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  2.9e+03
                        O2   3.62e+03
                        N2   1.76e+04
    [1] makeup_MEA
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): MEA  2.82
    [2] makeup_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water  22.3
    outs...
    [0] vent
        phase: 'l', T: 313.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  290
                        O2   3.62e+03
                        N2   1.76e+04
    [1] CO2
        phase: 'l', T: 313.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  2.61e+03
    >>> U1.results()
    Amine absorption                            Units       U1
    Power               Rate                       kW 1.23e+03
                        Cost                   USD/hr     96.4
    Low pressure steam  Duty                    kJ/hr 4.36e+08
                        Flow                  kmol/hr 1.12e+04
                        Cost                   USD/hr 2.67e+03
    Design              Total flow            kmol/hr 2.41e+04
                        CO2 flow              kmol/hr 2.61e+03
    Purchase cost       Makeup tank               USD  2.5e+05
                        Cooler                    USD 9.78e+04
                        Cross heat exchanger      USD 2.48e+06
                        Reboiler                  USD 5.76e+05
                        Condenser                 USD 2.94e+05
                        Pumps                     USD  1.2e+05
                        Stripper                  USD 4.35e+06
                        Absorber                  USD 5.23e+06
    Total purchase cost                           USD 1.34e+07
    Utility cost                               USD/hr 2.77e+03
    
    References
    ----------
    .. [1] Karimi et al., Capital Costs and Energy Considerations of Different
        Alternative Stripper Configurations for Post Combustion CO2 Capture.
        Chemical Engineering Research and Design 2011, 89 (8), 1229–1236.
        https://doi.org/10.1016/j.cherd.2011.03.005.
    
    .. [2] Carminati et al., Bioenergy and Full Carbon Dioxide Sinking in
        Sugarcane-Biorefinery with Post-Combustion Capture and Storage:
        Techno-Economic Feasibility. Applied Energy 2019, 254, 113633.
        https://doi.org/10.1016/j.apenergy.2019.113633.
        
    .. [3] Ramezan et al., Carbon Dioxide Capture from Existing Coal-Fired Power
        Plants; DOE/NETL-401/110907; National Energy Technology Laboratory, 2007.
    '''
    
    _N_ins = 3
    _N_outs = 2
    _N_heat_utilities = 1
    _units = {'Total flow': 'kmol/hr',
              'CO2 flow':'kmol/hr'}
    
    def __init__(self, ID='', ins=(), outs=(), thermo=None, *,
                 CO2_recovery=0.9, MEA_to_CO2=1.5, heat_ratio=3611):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.CO2_recovery = CO2_recovery
        self.MEA_to_CO2 = MEA_to_CO2
        self.heat_ratio = heat_ratio
    
    def _run(self):
        flue_gas, MEA, water = self.ins
        vent, CO2 = self.outs
        vent.copy_like(flue_gas)
        CO2.imol['CO2'] = self.CO2_recovery * flue_gas.imol['CO2']
        vent.imol['CO2'] = flue_gas.imol['CO2'] - CO2.imol['CO2']
        MEA.imass['MEA'] = self.MEA_to_CO2 * CO2.imass['CO2']/1000
        water.imass['Water'] = MEA.imass['MEA'] / 0.3 * (1-0.3)
        vent.T = CO2.T = 273.15 + 40
        
    def _design(self):
        self.design_results['Total flow'] = self.ins[0].F_mol
        self.design_results['CO2 flow'] = self.outs[1].F_mol
        duty = self.heat_ratio * self.outs[1].F_mass
        self.heat_utilities[0](duty, T_in=self.ins[0].T)
    

@cost(basis='Total flow', ID='Compressor', units='kmol/hr',
      # 17264.08 for (613+39.5) metric tonne/hr CO2
      kW=17264.08/((613+39.5)*1000/44)*(24123*0.1186),
      cost=5.08e6, S=24123, CE=CEPCI_by_year[2009], n=0.6, BM=3.3)
class CO2Compression(Unit):
    '''
    Create a CO2Compression unit for compression of gas CO2 into liquid CO2 [1]_.
    Cost is extrapolated via the 6/10th rule [2]_. Pump power usage is
    based on [2]_. Bare module factor is based on BioSTEAM pump.
    
    Parameters
    ----------
    ins : stream
        CO2 gas.
    outs : stream
        Compressed CO2.
    
    Examples
    --------
    
    >>> import biosteam as bst
    >>> import thermosteam as tmo
    >>> tmo.settings.set_thermo((tmo.Chemical('CO2'),))
    >>> gaseous_CO2 = tmo.Stream('gaseous_CO2',
    ...                          CO2=2.61e+03,
    ...                          units='kmol/hr')
    >>> U1 = bst.units.CO2Compression('U1', ins=gaseous_CO2, outs='compressed_CO2')
    >>> U1.simulate()
    >>> U1.show()
    CO2Compression: U1
    ins...
    [0] gaseous_CO2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  2.61e+03
    outs...
    [0] compressed_CO2
        phase: 'l', T: 303.15 K, P: 1.1e+07 Pa
        flow (kmol/hr): CO2  2.61e+03
    >>> U1.results()
    CO2Compression                    Units       U1
    Power               Rate             kW      360
                        Cost         USD/hr     28.2
    Design              Total flow  kmol/hr 2.61e+03
    Purchase cost       Compressor      USD 1.45e+06
    Total purchase cost                 USD 1.45e+06
    Utility cost                     USD/hr     28.2

    References
    ----------
    .. [1] Karimi et al., Capital Costs and Energy Considerations of Different
        Alternative Stripper Configurations for Post Combustion CO2 Capture.
        Chemical Engineering Research and Design 2011, 89 (8), 1229–1236.
        https://doi.org/10.1016/j.cherd.2011.03.005.
    
    .. [2] Carminati et al., Bioenergy and Full Carbon Dioxide Sinking in
        Sugarcane-Biorefinery with Post-Combustion Capture and Storage:
        Techno-Economic Feasibility. Applied Energy 2019, 254, 113633.
        https://doi.org/10.1016/j.apenergy.2019.113633.
    '''

    _N_ins = 1
    _N_outs = 1
    _units = {'Total flow':'kmol/hr'}
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        self.outs[0].T = 273.15 + 30
        self.outs[0].P = 110 * 1e5
        self.outs[0].phase = 'l'
    
    def _design(self):
        self.design_results['Total flow'] = self.outs[0].F_mol































    