# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

__all__ = ('SLLECentrifuge', 'SolidLiquidsSplitCentrifuge',)

from .design_tools import CEPCI_by_year
from thermosteam import separations as sep
from .decorators import cost, copy_algorithm
from . import Unit


@cost('Flow rate', S=1725.61, units='L/min',
    CE=CEPCI_by_year[2007], cost=849000., n=0.6, kW=0.07, ub=2000., BM=2.03,
    N='Number of centrifuges')
class SLLECentrifuge(Unit):
    """
    Create a SLLECentrifuge object that separates the feed into solid, oil, and
    aqueous phases.
    
    Parameters
    ----------
    ins : stream
        feed
    outs : stream sequence
        [0] Oil fluid.
        [1] Aqueous fluid.
        [2] Solids.
    solids_split : dict[str, float]
        Splits to 2nd outlet stream.
    top_chemical : str, optional
        Identifier of chemical that will be favored in the oil phase.
        If none given, the oil phase will the lightest and the aqueous
        phase will be the heaviest.
    efficiency : float, optional
        Fraction of feed in liquid-liquid equilibrium.
        The rest of the feed is divided equally between phases.
        Defaults to 1.0.
    moisture_content : float, optional
        Moisture content of solids. Defaults to 0.5.
    
    Notes
    -----
    Cost algorithm is based on a 3-phase decanter centrifuge from 
    a conventional dry-grind corn ethanol plant that separates
    aqueous, oil, and solid fractions (i.e. DDGS) from the bottoms product 
    of the beer column [1]_.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Hexane', bst.Chemical('Solids', search_db=False, default=True, phase='s')], cache=True)
    >>> feed = bst.Stream('feed', Water=100, Hexane=100, Solids=10)
    >>> C1 = bst.SLLECentrifuge('C1', feed, ['oil', 'aqueous', 'solids'], solids_split={'Solids':1.0})
    >>> C1.simulate()
    >>> C1.show()
    SLLECentrifuge: C1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   100
                        Hexane  100
                        Solids  10
    outs...
    [0] oil
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   0.791
                        Hexane  100
    [1] aqueous
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   98.7
                        Hexane  0.015
    [2] solids
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water   0.555
                        Solids  10
    
    >>> C1.results()
    3-Phase decanter centrifuge                       Units       C1
    Power               Rate                             kW   0.0101
                        Cost                         USD/hr 0.000792
    Design              Flow rate                     L/min      250
                        Number of centrifuges                      1
    Purchase cost       3-Phase decanter centrifuge     USD 2.88e+05
    Total purchase cost                                 USD 2.88e+05
    Utility cost                                     USD/hr 0.000792
    
    References
    ----------
    .. [1] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
        Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
        Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
        https://doi.org/10.1016/j.indcrop.2005.08.004.

    """
    line = '3-Phase decanter centrifuge'
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 0
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 solids_split, top_chemical=None, efficiency=1.0,
                 moisture_content=0.5):
        Unit.__init__(self, ID, ins, outs, thermo)
        
        # [ChemicalIndexer] Splits to 0th outlet stream.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        #: [str] Identifier of chemical that will be favored in the "liquid" phase.
        #: If none given, the "liquid" phase will the lightest and the "LIQUID"
        #: phase will be the heaviest.
        self.top_chemical = top_chemical
        
        #: Fraction of feed in liquid-liquid equilibrium.
        #: The rest of the feed is divided equally between phases.
        self.efficiency = efficiency 
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split, only moisture content'

    def _run(self):
        top, bottom, solids = self.outs
        sep.lle(self.ins[0], top, bottom, self.top_chemical, self.efficiency)
        sep.split(bottom, solids, bottom, self.solids_split)
        sep.adjust_moisture_content(solids, bottom, self.moisture_content)


@copy_algorithm(SLLECentrifuge, run=False)
class SolidLiquidsSplitCentrifuge(Unit):
    """
    Create a SolidLiquidsSplitCentrifuge object that separates the feed into solid, oil, and
    aqueous phases.
    
    Parameters
    ----------
    ins : stream
        feed
    outs : stream sequence
        [0] Oil fluid.
        [1] Aqueous fluid.
        [2] Solids.
    aqueous_split : dict[str, float]
        Splits to [0] outlet stream.
    solids_split : dict[str, float]
        Splits to [2] outlet stream.
    moisture_content : float, optional
        Moisture content of solids. Defaults to 0.5.
    
    Notes
    -----
    Cost algorithm is based on a 3-phase decanter centrifuge from 
    a conventional dry-grind corn ethanol plant that separates
    aqueous, oil, and solid fractions (i.e. DDGS) from the bottoms product 
    of the beer column [1]_.
    
    The unit operation first splits the feed to the aqueous and oil fluids, 
    then fractionates the solids from the aqueous phase.
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Hexane', bst.Chemical('Solids', search_db=False, default=True, phase='s')], cache=True)
    >>> feed = bst.Stream('feed', Water=100, Hexane=100, Solids=10)
    >>> C1 = bst.SolidLiquidsSplitCentrifuge(
    ...     'C1', feed, ['oil', 'aqueous', 'solids'],
    ...     solids_split={'Solids':1.0},
    ...     aqueous_split={'Water':0.99, 'Hexane':0.01, 'Solids': 0.9}
    ... )
    >>> C1.simulate()
    >>> C1.show(flow='kg/hr')
    SolidLiquidsSplitCentrifuge: C1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   1.8e+03
                      Hexane  8.62e+03
                      Solids  10
    outs...
    [0] oil
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   18
                      Hexane  8.53e+03
                      Solids  1
    [1] aqueous
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   1.77e+03
                      Hexane  86.2
    [2] solids
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): Water   9
                      Solids  9
    
    >>> C1.results()
    3-Phase decanter centrifuge                       Units       C1
    Power               Rate                             kW   0.0101
                        Cost                         USD/hr 0.000792
    Design              Flow rate                     L/min      250
                        Number of centrifuges                      1
    Purchase cost       3-Phase decanter centrifuge     USD 2.88e+05
    Total purchase cost                                 USD 2.88e+05
    Utility cost                                     USD/hr 0.000792
    
    References
    ----------
    .. [1] Kwiatkowski, J. R.; McAloon, A. J.; Taylor, F.; Johnston, D. B. 
        Modeling the Process and Costs of Fuel Ethanol Production by the Corn 
        Dry-Grind Process. Industrial Crops and Products 2006, 23 (3), 288–296.
        https://doi.org/10.1016/j.indcrop.2005.08.004.
    
    """
    line = SLLECentrifuge.line
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 0
    
    @property
    def solids_split(self):
        return self._solids_isplit.data
    @property
    def solids_isplit(self):
        return self._solids_isplit
    
    @property
    def aqueous_split(self):
        return self._aqueous_isplit.data
    @property
    def aqueous_isplit(self):
        return self._aqueous_isplit
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 aqueous_split, solids_split, moisture_content=0.5):
        Unit.__init__(self, ID, ins, outs, thermo)
        
        # [ChemicalIndexer] Splits to 1st outlet stream (aqueous/heavy phase).
        self._aqueous_isplit = self.thermo.chemicals.isplit(aqueous_split)
        
        # [ChemicalIndexer] Splits to 2th outlet stream (solids) from aqueous/heavy phase.
        self._solids_isplit = self.thermo.chemicals.isplit(solids_split)
        
        #: Moisture content of retentate
        self.moisture_content = moisture_content
        assert self._solids_isplit['7732-18-5'] == 0, 'cannot define water split to solids, only moisture content'

    def _run(self):
        oil, aqueous, solids = self.outs
        sep.split(*self.ins, aqueous, oil, self.aqueous_split)
        sep.split(aqueous, solids, aqueous, self.solids_split)
        sep.adjust_moisture_content(solids, aqueous, self.moisture_content)
        
        
        