# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import Facility
from ..decorators import cost

__all__ = ('ProcessWaterCenter',)

@cost('Makeup water flow rate', 'Makeup water pump',
      CE=551, kW=20*0.7457, cost=6864, S=155564, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Process water pump',
      CE=551, kW=75*0.7457, cost=15292, S=518924, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Tank',
      CE=522, cost=250e3, S=451555, n=0.7, BM=1.7)
class ProcessWaterCenter(Facility):
    """
    Create a ProcessWaterCenter object that takes care of balancing the amount
    of make-up water required for the process. The capital cost and power
    are based on the flow rate of process and make-up water as in [1]_.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Clean water.
        
        [1] Make-up water.
        
        [2] Recycled process water
    outs : stream
        [0] Process water.
        
        [1] Waste.
    makeup_water_streams : streams, optional
        All fresh makeup water streams (must be a subset of `process_water_streams`).
    process_water_streams : streams, optional
        All process water streams (including makeup water streams).
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    
    """
    ticket_name = 'PWC'
    network_priority = 2
    _N_ins = 3
    _N_outs = 2
    _units = {'Makeup water flow rate': 'kg/hr',
              'Process water flow rate': 'kg/hr'}
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 makeup_water_streams=None,
                 process_water_streams=None):
        Facility.__init__(self, ID, ins, outs, thermo)
        self.makeup_water_streams = makeup_water_streams
        self.process_water_streams = process_water_streams
    
    def _assert_compatible_property_package(self): pass
    
    def update_process_water(self):
        process_water_streams = self.process_water_streams
        s_process, _ = self.outs
        process_water = sum([stream.imol['7732-18-5'] 
                             for stream in process_water_streams])
        s_process.imol['7732-18-5'] = process_water

    def update_makeup_water(self):
        makeup_water_streams = self.makeup_water_streams
        s_recycle, s_makeup, _ = self.ins
        water = sum([stream.imol['7732-18-5'] for stream in makeup_water_streams]) - s_recycle.imol['7732-18-5']
        s_makeup.imol['7732-18-5'] = max(water, 0)

    def _run(self): 
        self.update_process_water()
        self.update_makeup_water()
        s_recycle, s_makeup, s_recycle_process = self._ins
        s_process, s_waste = self.outs
        makeup_water = s_makeup.imol['7732-18-5']
        recycle_water = sum([i.imol['7732-18-5'] for i in (s_recycle, s_recycle_process) if i])
        process_water = s_process.imol['7732-18-5']
        waste_water = recycle_water + makeup_water - process_water
        if waste_water < 0.:
            s_makeup.imol['7732-18-5'] -= waste_water
            waste_water = 0.
        s_waste.imol['7732-18-5'] = waste_water
        Design = self.design_results
        Design['Process water flow rate'] = (process_water + waste_water) * 18.015
        Design['Makeup water flow rate'] = makeup_water * 18.015
        
        