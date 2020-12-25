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
from ._liquids_mixing_tank import LiquidsMixingTank
from ._liquids_settler import LiquidsSettler, LLESettler, LiquidsSplitSettler, LiquidsPartitionSettler

__all__ = ('MixerSettler',
           'MultiStageMixerSettlers')

class MixerSettler(bst.Unit):
    """
    Create a MixerSettler object that models liquid-liquid extraction using
    a mixing tank and a settler tank.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] feed.
        * [1] solvent.
    outs : stream sequence
        * [0] raffinate
        * [1] extract
    carrier_chemical : str, optional
        Name of main chemical in the feed (which is not selectively extracted by the solvent).
        Defaults to chemical with highest molar fraction in the feed.
    mixer_data : dict, optional
        Arguments to initialize the "mixer" attribute, a :class:`~biosteam.units.LiquidsMixingTank` object.
    settler_data : dict, optional
        Arguments to initialize the "settler" attribute, a :class:`~biosteam.units.LiquidsSettler` object.
    
    Examples
    --------
    Simulate by rigorous LLE:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MS1 = bst.MixerSettler('MS1', ins=(feed, solvent), outs=('raffinate', 'extract'))
    >>> MS1.simulate()
    >>> MS1.show()
    MixerSettler: MS1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     500
                        Methanol  50
    [1] solvent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octanol  500
    outs...
    [0] raffinate
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     414
                        Methanol  16.7
                        Octanol   0.117
    [1] extract
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     86.3
                        Methanol  33.3
                        Octanol   500
    >>> MS1.results()
    Mixer settler                                  Units         MS1
    Power               Rate                          kW        1.98
                        Cost                      USD/hr       0.155
    Design              Mixer - Volume               m^3        1.98
                        Mixer - Power                 hp        2.66
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        1.36
                        Mixer - Diameter              ft        1.36
                        Mixer - Weight                lb        91.4
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length              ft        12.6
                        Settler - Diameter            ft        3.15
                        Settler - Weight              lb    1.44e+03
                        Settler - Wall thickness      in        0.25
    Purchase cost       Mixer                        USD    1.16e+04
                        Settler                      USD    2.89e+03
    Total purchase cost                              USD    1.45e+04
    Utility cost                                  USD/hr       0.155
    
    Simulate with user defined partition coefficients:
    
    >>> import biosteam as bst
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MS1 = bst.MixerSettler('MS1', 
    ...    ins=(feed, solvent), outs=('raffinate', 'extract'),
    ...    model='partition coefficients',
    ...    settler_data={
    ...        'partition_coefficients': np.array([6.894, 0.7244, 3.381e-04]),
    ...        'partion_IDs': ('Water', 'Methanol', 'Octanol'),
    ...    },
    ... )
    >>> MS1.simulate()
    >>> MS1.show()
    MixerSettler: MS1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     500
                        Methanol  50
    [1] solvent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octanol  500
    outs...
    [0] raffinate
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     414
                        Methanol  16.7
                        Octanol   0.117
    [1] extract
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     86.3
                        Methanol  33.3
                        Octanol   500
    >>> MS1.results()
    Mixer settler                                  Units         MS1
    Power               Rate                          kW        1.98
                        Cost                      USD/hr       0.155
    Design              Mixer - Volume               m^3        1.98
                        Mixer - Power                 hp        2.66
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        1.36
                        Mixer - Diameter              ft        1.36
                        Mixer - Weight                lb        91.4
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length              ft        12.6
                        Settler - Diameter            ft        3.15
                        Settler - Weight              lb    1.44e+03
                        Settler - Wall thickness      in        0.25
    Purchase cost       Mixer                        USD    1.16e+04
                        Settler                      USD    2.89e+03
    Total purchase cost                              USD    1.45e+04
    Utility cost                                  USD/hr       0.155
    
    """
    _N_ins = 2
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('mixer', 'settler')
    _graphics = mixer_settler_graphics
    _units = {}
    for i,j in LiquidsMixingTank._units.items(): 
        _units['Mixer - ' + i] = j
    for i,j in LiquidsSettler._units.items(): 
        _units['Settler - ' + i] = j
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 carrier_chemical=None, mixer_data={}, settler_data={}, model="LLE"):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        
        #: [LiquidsMixingTank] Mixer portion of the mixer-settler.
        #: All data and settings for the design of the mixing tank are stored here.
        self.mixer = mixer = LiquidsMixingTank(None, None, (None,),
                                                   self.thermo, **mixer_data)
        self.multi_stream = multi_stream = mixer-0
        mixer._ins = self._ins
        model = model.lower()
        if model == 'lle':
            Settler = LLESettler
        elif model == 'split':
            Settler = LiquidsSplitSettler
        elif model == 'partition coefficients':
            Settler = LiquidsPartitionSettler
            
        #: [LiquidsSettler] Settler portion of the mixer-settler.
        #: All data and settings for the design of the settler are stored here.
        self.settler = Settler(None, multi_stream, None, self.thermo, **settler_data)
        self.settler._outs = self._outs
        self.power_utility = mixer.power_utility
        
        #: [str] ID of carrier component in the feed.
        self.carrier_chemical = carrier_chemical 
    
    @property
    def feed(self):
        """[Stream] Feed with solute being extracted and carrier."""
        return self._ins[0]
    @feed.setter
    def feed(self, stream):
        self._ins[0] = stream

    @property
    def solvent(self):
        """[Stream] Solvent to extract solute."""
        return self._ins[1]
    @solvent.setter
    def solvent(self, stream):
        self._ins[1] = stream
        
    @property
    def raffinate(self):
        """[Stream] Raffinate after extraction."""
        return self._outs[0]
    @raffinate.setter
    def raffinate(self, stream):
        self._outs[0] = stream
        
    @property
    def extract(self):
        """[Stream] Extract with solvent."""
        return self._outs[1]
    @extract.setter
    def extract(self, stream):
        self._outs[1] = stream

    def _run(self):
        self.mixer._run()
        self.settler.top_chemical = self.carrier_chemical or self.feed.main_chemical
        self.settler._run()
        
    def _design(self):
        mixer = self.mixer
        mixer._design()
        settler = self.settler
        settler._design()
        design_results = self.design_results
        for i,j in mixer.design_results.items():
            design_results['Mixer - ' + i] = j
        for i,j in settler.design_results.items():
            design_results['Settler - ' + i] = j
        
    def _cost(self):
        self.mixer._cost()
        self.settler._cost()
        
class MultiStageMixerSettlers(bst.Unit):
    """
    Create a MultiStageMixerSettlers object that models a counter-current system
    of mixer-settlers for liquid-liquid extraction.
    
    Parameters
    ----------
    ins : stream sequence
        * [0] feed.
        * [1] solvent.
    outs : stream sequence
        * [0] raffinate
        * [1] extract
    N_stages : int
        Number of stages.
    partition_data : {'IDs': tuple[str], 'K': 1d array}, optional
        IDs of chemicals in equilibrium and partition coefficients (molar 
        composition ratio of the raffinate over the extract). If given,
        The mixer-settlers will be modeled with these constants. Otherwise,
        partition coefficients are computed based on temperature and composition.
    carrier_chemical : str
        Name of main chemical in the feed (which is not selectively extracted by the solvent).
    mixer_data : dict
        Arguments to initialize the "mixer" attribute, a :class:`~biosteam.units.LiquidsMixingTank` object.
    settler_data : dict
        Arguments to initialize the "settler" attribute, a :class:`~biosteam.units.LiquidsSettler` object.
    
    Notes
    -----
    All mixer settlers are sized equally based on the assumption that the 
    volumetric flow rate of each phase does not change significantly across
    stages. 
    
    Examples
    --------
    Simulate by rigorous LLE:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), outs=('raffinate', 'extract'), N_stages=2)
    >>> MSMS1.simulate()
    >>> MSMS1.show()
    MultiStageMixerSettlers: MSMS1
    ins...
    [0] feed
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     500
                        Methanol  50
    [1] solvent
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Octanol  500
    outs...
    [0] raffinate
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     413
                        Methanol  8.4
                        Octanol   0.103
    [1] extract
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water     87.2
                        Methanol  41.6
                        Octanol   500
    >>> MSMS1.results()
    Multi stage mixer settlers                     Units       MSMS1
    Power               Rate                          kW        3.97
                        Cost                      USD/hr        0.31
    Design              Mixer - Volume               m^3        1.98
                        Mixer - Power                 hp        2.66
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        1.36
                        Mixer - Diameter              ft        1.36
                        Mixer - Weight                lb        91.4
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length              ft        12.6
                        Settler - Diameter            ft        3.15
                        Settler - Weight              lb    1.44e+03
                        Settler - Wall thickness      in        0.25
    Purchase cost       Mixers and agitators         USD    2.31e+04
                        Settlers                     USD    5.78e+03
    Total purchase cost                              USD    2.89e+04
    Utility cost                                  USD/hr        0.31
    
    Simulate with user defined partition coefficients:
    
    >>> import biosteam as bst
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
    >>> MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins=(feed, solvent), outs=('raffinate', 'extract'), N_stages=10,
    ...     partition_data={
    ...         'K': np.array([6.894, 0.7244, 3.381e-04]),
    ...         'IDs': ('Water', 'Methanol', 'Octanol'),
    ...         'phi': 0.4100271108219455 # Initial phase fraction guess. This is optional.
    ...     }
    ... )
    >>> MSMS1.simulate()
    >>> MSMS1.raffinate.show()
    Stream: raffinate from <MultiStageMixerSettlers: MSMS1>
     phase: 'l', T: 298.15 K, P: 101325 Pa
     flow (kmol/hr): Water     4.13e+03
                     Methanol  1.27
                     Octanol   1.51
    
    >>> MSMS1.results()
    Multi stage mixer settlers                     Units       MSMS1
    Power               Rate                          kW         198
                        Cost                      USD/hr        15.5
    Design              Mixer - Volume               m^3        19.8
                        Mixer - Power                 hp        26.6
                        Mixer - Vessel type                 Vertical
                        Mixer - Length                ft        2.93
                        Mixer - Diameter              ft        2.93
                        Mixer - Weight                lb         424
                        Mixer - Wall thickness        in        0.25
                        Settler - Vessel type             Horizontal
                        Settler - Length              ft        39.9
                        Settler - Diameter            ft        9.96
                        Settler - Weight              lb    2.53e+04
                        Settler - Wall thickness      in       0.438
    Purchase cost       Mixers and agitators         USD    3.27e+05
                        Settlers                     USD    3.64e+04
    Total purchase cost                              USD    3.64e+05
    Utility cost                                  USD/hr        15.5
    
    """
    _units = MixerSettler._units
    _N_ins = 2
    _N_outs = 2
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *, N_stages,
                 partition_data=None, carrier_chemical=None,
                 mixer_data={}, settler_data={}):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        
        #: [str] Name of main chemical in the feed (which is not extracted by 
        #: the solvent).
        self.carrier_chemical = carrier_chemical
        
        #: [int] Number of stages.
        self.N_stages = N_stages
        
        # {'IDs': tuple[str], 'K': 1d array}
        # IDs of chemicals in equilibrium and partition coefficients. If given,
        # The mixer-settlers will be modeled with these constants. Otherwise,
        # partition coefficients are computed based on temperature and composition.
        self.partition_data = partition_data
        
        #: [LiquidsMixingTank] Used to design all mixing tanks. 
        #: All data and settings for the design of mixing tanks are stored here.
        self.mixer = mixer = LiquidsMixingTank(None, None, (None,),
                                                   self.thermo, **mixer_data)
        mixer._ins = self._ins
        
        #: [LiquidsSettler] Used to design all settlers.
        #: All data and settings for the design of settlers are stored here.
        self.settler = LiquidsSettler(None, mixer-0, None, self.thermo, **settler_data)
        
        self.settler._outs = self._outs
        self.reset_cache()
    
    feed = MixerSettler.feed
    solvent = MixerSettler.solvent
    raffinate = MixerSettler.raffinate
    extract = MixerSettler.extract
    
    def _setup(self):
        args = (self.stages, self.feed, self.solvent, self.carrier_chemical)
        if args != self._last_args:
            self.stages = tmo.separations.MultiStageLLE(
                self.N_stages, self.feed, self.solvent, self.carrier_chemical, 
                self.thermo, self.partition_data
            )
            self._last_args = args
    
    def reset_cache(self):
        self.stages = None
        self._last_args = None
        
    def _run(self):
        stages = self.stages
        stages.simulate_multi_stage_lle_without_side_draws()
        self.extract.copy_flow(stages[0].extract)
        self.raffinate.copy_flow(stages[-1].raffinate)
        
    def _design(self):
        mixer = self.mixer
        mixer._run()
        mixer._design()
        settler = self.settler
        settler._design()
        design_results = self.design_results
        for i,j in mixer.design_results.items():
            design_results['Mixer - ' + i] = j
        for i,j in settler.design_results.items():
            design_results['Settler - ' + i] = j
        
    def _cost(self):
        purchase_costs = self.purchase_costs
        N_stages = self.N_stages
        mixer = self.mixer
        settler = self.settler
        mixer._cost()
        settler._cost()
        self.power_utility.copy_like(mixer.power_utility)
        self.power_utility.scale(N_stages)
        purchase_costs['Mixers and agitators'] = N_stages * mixer.purchase_cost
        purchase_costs['Settlers'] = N_stages * settler.purchase_cost
        
    @property
    def installed_cost(self):
        N_stages = self.N_stages
        return N_stages * (self.mixer.installed_cost + self.settler.installed_cost)
        