# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
cost = bst.decorators.cost

__all__ = ('ProcessWaterCenter',)

@cost('Makeup water flow rate', 'Makeup water pump',
      CE=551, kW=20*0.7457, cost=6864, S=155564, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Process water pump',
      CE=551, kW=75*0.7457, cost=15292, S=518924, n=0.8, BM=3.1)
@cost('Process water flow rate', 'Tank',
      CE=522, cost=250e3, S=451555, n=0.7, BM=1.7)
class ProcessWaterCenter(bst.Facility):
    """
    Create a ProcessWaterCenter object that takes care of balancing the amount
    of make-up process and reverse osmosis (RO) -grade water required for the
    process. The capital cost and power are based on the flow rate of process
    and make-up water as given in[1]_.
    
    Parameters
    ----------
    ins : 
        [0] Recycled RO-grade water.
        
        [1] Make-up RO-grade water.
        
        [2] Recycled process water.  
        
        [3] Make-up process water.
    outs : 
        [0] RO-grade water.
        
        [1] Process water.
        
        [2] Excess water.
    makeup_water_streams : List[Stream], optional
        All inlet RO-grade water streams.
        Defaults to boiler and cooling tower make-up water streams at run time.
    process_water_streams : List[Stream], optional
        All inlet process water streams (excluding makeup water streams).
        Defaults to all fresh process water streams within the system at 
        run time.
    reverse_osmosis_water_price : float, optional
        Defaults to 5.6e-4 USD/kg.
    process_water_price : float, optional
        Defaults to 2.7e-4 USD/kg.
        
    Notes
    -----
    Default prices for the RO-grade and process water are 0.56 and 0.27 USD/m3 
    as given in Table 17.1 of [2]_.
    
    References
    ----------
    .. [1] Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., Aden, A.,
        Dudgeon, D. (2011). Process Design and Economics for Biochemical 
        Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid 
        Pretreatment and Enzymatic Hydrolysis of Corn Stover
        (No. NREL/TP-5100-47764, 1013269). https://doi.org/10.2172/1013269
    .. [2] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
        R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Cost Accounting and Capital Cost Estimation (Chapter 16)
    
    """
    ticket_name = 'PWC'
    network_priority = 2
    _N_ins = 4
    _N_outs = 3
    _units = {'Makeup water flow rate': 'kg/hr',
              'Process water flow rate': 'kg/hr'}
    def _init(self, 
            makeup_water_streams=None,
            process_water_streams=None,
            reverse_osmosis_water_price=None,
            process_water_price=None
        ):
        if process_water_streams and makeup_water_streams:
            process_water_streams = list(process_water_streams)
            for i in makeup_water_streams:
                if i in process_water_streams: process_water_streams.remove(i)
        self.makeup_water_streams = makeup_water_streams
        self.process_water_streams = process_water_streams
        self.define_utility('Reverse osmosis water', self.makeup_reverse_osmosis_grade_water)
        self.define_utility('Process water', self.makeup_process_water)
        if reverse_osmosis_water_price: self.reverse_osmosis_water_price = reverse_osmosis_water_price
        if process_water_price: self.process_water_price = process_water_price
    
    @property
    def reverse_osmosis_water_price(self):
        """[Float] Price of reverse osmosis-grade water, same as `bst.stream_utility_prices['Reverse osmosis water']`."""
        return bst.stream_utility_prices['Reverse osmosis water']
    
    @reverse_osmosis_water_price.setter
    def reverse_osmosis_water_price(self, new_price):
        bst.stream_utility_prices['Reverse osmosis water'] = new_price
    
    @property
    def process_water_price(self):
        """[Float] Price of process water, same as `bst.stream_utility_prices['Process water']`."""
        return bst.stream_utility_prices['Process water']
    
    @process_water_price.setter
    def process_water_price(self, new_price):
        bst.stream_utility_prices['Process water'] = new_price
    
    @property
    def recycled_reverse_osmosis_grade_water(self):
        return self.ins[0]
    
    @property
    def makeup_reverse_osmosis_grade_water(self):
        return self.ins[1]
    
    @property
    def recycled_process_water(self):
        return self.ins[2]
    
    @property
    def makeup_process_water(self):
        return self.ins[3]
        
    @property
    def reverse_osmosis_grade_water(self):
        return self.outs[0]
    
    @property
    def process_water(self):
        return self.outs[1]
    
    @property
    def excess_water(self):
        return self.outs[2]
    
    def _assert_compatible_property_package(self): pass
    
    def update_process_water(self):
        process_water_streams = self.process_water_streams
        makeup_water_streams = set(self.makeup_water_streams)
        if process_water_streams is None:
            self.process_water_streams = process_water_streams = [i for i in bst.get_fresh_process_water_streams() if i not in makeup_water_streams]
        process_water = sum([i.imol['7732-18-5'] for i in process_water_streams])
        self.process_water.imol['7732-18-5'] = process_water

    def update_reverse_osmosis_grade_water(self):
        makeup_water_streams = self.makeup_water_streams
        if makeup_water_streams is None: 
            self.makeup_water_streams = makeup_water_streams = [
                i.makeup_water for i in self.system.facilities
                if hasattr(i, 'makeup_water')
            ]
        self.reverse_osmosis_grade_water.imol['7732-18-5'] = sum([stream.imol['7732-18-5'] for stream in makeup_water_streams])

    def _run(self): 
        self.update_reverse_osmosis_grade_water()
        self.update_process_water()
        s_recycle_rev_os, s_rev_os_makeup, s_recycle_process, s_process_makeup = self._ins
        rrogw = s_recycle_rev_os.imol['7732-18-5']
        rpw = s_recycle_process.imol['7732-18-5']
        pw = self.process_water.imol['7732-18-5']
        rogw = self.reverse_osmosis_grade_water.imol['7732-18-5']
        mrogw = rogw - rrogw
        if mrogw >= 0.:
            s_rev_os_makeup.imol['7732-18-5'] = mrogw
        else:
            s_rev_os_makeup.imol['7732-18-5'] = 0.
            rpw -= mrogw # Unused RO-grade water is added to recycled process water
        mpw = pw - rpw
        if mpw >= 0.:
            s_process_makeup.imol['7732-18-5'] = mpw
            self.excess_water.imol['7732-18-5'] = 0.
        else:
            s_process_makeup.imol['7732-18-5'] = 0.
            self.excess_water.imol['7732-18-5'] = -mpw # Excess recycled water
        Design = self.design_results
        Design['Process water flow rate'] = (
            self.excess_water.imol['7732-18-5'] + self.process_water.imol['7732-18-5']
        ) * 18.015
        Design['Makeup water flow rate'] = (
            s_process_makeup.imol['7732-18-5']
        ) * 18.015
        
        