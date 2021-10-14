# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .splitting import Splitter
from .decorators import cost

__all__ = ('MolecularSieve',)

@cost('Flow rate', 'Column', kW=151, BM=1.8,
      cost=2601000, CE=521.9, S=22687, n=0.6)
class MolecularSieve(Splitter):
    """
    Create an ethanol/water molecular sieve for bioethanol plants.
    The molecular sieve is modeled as a component wise separator. Costing
    is based on scaling by the 6/10ths rule from an NREL TEA report [1]_.
    
    Parameters
    ----------
    ins : stream
        * [0] Feed (gas)
    outs : stream sequence
        * [0] Split stream (gas)
        * [1] Remainder stream (gas)
    split : array_like
            Componentwise split to the 0th output stream
    
    Examples
    --------
    >>> from biosteam import Stream, settings
    >>> from biosteam.units import MolecularSieve
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = Stream('feed', flow=(75.7, 286), T=351.39, phase='g')
    >>> bp = feed.bubble_point_at_T()
    >>> feed.T = bp.T
    >>> MS1 = MolecularSieve('MS1', ins=feed,
    ...                      outs=('ethanol_rich', 'water_rich'),
    ...                      split=dict(Water=0.160,
    ...                                 Ethanol=0.925))
    >>> MS1.simulate()
    >>> MS1.show(T='degC', P='atm', composition= True)
    MolecularSieve: MS1
    ins...
    [0] feed
        phase: 'g', T: 78.24 degC, P: 1 atm
        composition: Water    0.209
                     Ethanol  0.791
                     -------  362 kmol/hr
    outs...
    [0] ethanol_rich
        phase: 'g', T: 78.24 degC, P: 1 atm
        composition: Water    0.0438
                     Ethanol  0.956
                     -------  277 kmol/hr
    [1] water_rich
        phase: 'g', T: 78.24 degC, P: 1 atm
        composition: Water    0.748
                     Ethanol  0.252
                     -------  85 kmol/hr
    
    >>> MS1.results()
    Molecular sieve                  Units       MS1
    Power               Rate            kW      14.2
                        Cost        USD/hr      1.11
    Low pressure steam  Duty         kJ/hr  3.21e+06
                        Flow       kmol/hr      82.7
                        Cost        USD/hr      19.7
    Cooling water       Duty         kJ/hr -1.18e+05
                        Flow       kmol/hr      80.9
                        Cost        USD/hr    0.0395
    Design              Flow rate    kg/hr  2.13e+03
    Purchase cost       Column         USD  6.85e+05
    Total purchase cost                USD  6.85e+05
    Utility cost                    USD/hr      20.8
    
    References
    ----------
    .. [1] Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and
        Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L.
        Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy
        Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof,
        M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle,
        Washington and Atlanta, Georgia)
    
    """
    _N_heat_utilities = 2
    _units = {'Flow rate': 'kg/hr'}
    def __init__(self, ID='', ins=None, outs=(), *, order=None, split,
                 P=101325, approx_duty=True):
        Splitter.__init__(self, ID, ins, outs, order=order, split=split)
        self.P = 101325
        self.approx_duty = approx_duty
        
    def _run(self):
        Splitter._run(self)
        for i in self.outs: i.P = self.P

    def _design(self):
        self.design_results['Flow rate'] = flow = self._outs[1].F_mass
        if self.approx_duty:
            T = self.outs[0].T
            hu1, hu2 = self.heat_utilities
            hu1(1429.65*flow, T)
            hu2(-55.51*flow, T)

