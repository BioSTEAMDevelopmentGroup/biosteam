# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
import thermosteam as tmo
from .._graphics import mixer_settler_graphics

__all__ = ('MixerSettler',
           'MultiStageMixerSettlers')

class MixerSettler(bst.Unit):
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('mixer', 'settler')
    _graphics = mixer_settler_graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 mixer_data={}, settler_data={}, model="LLE"):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.mixer = mixer = bst.LiquidsMixingTank(None, None, (None,),
                                                   self.thermo, **mixer_data)
        self.multi_stream = multi_stream = mixer-0
        mixer._ins = self._ins
        model = model.lower()
        if model == 'lle':
            Settler = bst.LLESettler
        elif model == 'split':
            Settler = bst.LiquidsSplitSettler
        elif model == 'partition coefficients':
            raise NotImplementedError("partition coefficient model not yet implemented in BioSTEAM")
            Settler = bst.LiquidsPartitionSettler
        self.settler = Settler(None, multi_stream, None, self.thermo, **settler_data)
        self.settler._outs = self._outs
        self.power_utility = mixer.power_utility
    
    @property
    def feed(self):
        return self._ins[0]
    @feed.setter
    def feed(self, stream):
        self._ins[0] = stream

    @property
    def solvent(self):
        return self._ins[1]
    @solvent.setter
    def solvent(self, stream):
        self._ins[1] = stream
        
    @property
    def raffinate(self):
        return self._outs[0]
    @raffinate.setter
    def raffinate(self, stream):
        self._outs[0] = stream
        
    @property
    def extract(self):
        return self._outs[1]
    @extract.setter
    def extract(self, stream):
        self._outs[1] = stream

    def _run(self):
        self.mixer._run()
        self.settler._run()
        
    def _design(self):
        self.mixer._design()
        self.settler._design()
        
    def _cost(self):
        self.mixer._cost()
        self.settler._cost()
        
class MultiStageMixerSettlers(bst.Unit):
    """
    This unit is not yet finished.
    
    Examples
    --------
    Simulate by rigorous LLE:
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream(Water=100, Methanol=10)
    >>> solvent = bst.Stream(Octanol=100)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), N_stages=2)
    >>> MSMS1.simulate()
    >>> MSMS1.show()
    MultiStageMixerSettlers: MSMS1
    ins...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     100
                        Methanol  10
    [1] s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octanol  100
    outs...
    [0] s3
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     82.6
                        Methanol  1.68
                        Octanol   0.0206
    [1] s4
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     17.4
                        Methanol  8.32
                        Octanol   100
    
    Simulate with user defined partition coefficients:
    >>> import biosteam as bst
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream(Water=100, Methanol=10)
    >>> solvent = bst.Stream(Octanol=100)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), N_stages=10,
    ...     partition_data={
    ...         'K': np.array([6.894, 0.7244, 3.381e-04]),
    ...         'IDs': ('Water', 'Methanol', 'Octanol'),
    ...         'phi': 0.4100271108219455
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> MSMS1.show()
    MultiStageMixerSettlers: MSMS1
    ins...
    [0] s1
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     100
                        Methanol  10
    [1] s2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octanol  100
    outs...
    [0] s3
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     82.6
                        Methanol  0.0295
                        Octanol   0.0234
    [1] s4
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     17.4
                        Methanol  9.97
                        Octanol   100
    """
    _N_ins = 2
    _N_outs = 2
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *, N_stages,
                 partition_data=None, carrier_chemical=None):
        super().__init__(ID, ins, outs, thermo)
        self.carrier_chemical = carrier_chemical
        self.N_stages = N_stages
        self.partition_data = partition_data
    
    feed = MixerSettler.feed
    solvent = MixerSettler.solvent
    raffinate = MixerSettler.raffinate
    extract = MixerSettler.extract
    
    def _setup(self):
        if not hasattr(self, 'stages'):
            self.stages = tmo.separations.MultiStageLLE(
                self.N_stages, self.feed, self.solvent, self.carrier_chemical, 
                self.thermo, self.partition_data
            )
    
    def reset_cache(self):
        del self.stages
        
    def _run(self):
        stages = self.stages
        stages.simulate_multi_stage_lle_without_side_draws()
        self.extract.copy_flow(stages[0].extract)
        self.raffinate.copy_flow(stages[-1].raffinate)
        