# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Created on Sat Aug 22 21:58:19 2020
@author: sarangbhagwat and yoelcp
"""
from .. import Facility
import biosteam as bst
from .hxn_synthesis import synthesize_network, StreamLifeCycle

__all__ = ('HeatExchangerNetwork',)


class HeatExchangerNetwork(Facility):
    """
    Create a HeatExchangerNetwork object that will perform a pinch analysis
    on the entire system's heating and cooling utility objects.
    
    Parameters
    ----------
    ID : str
        Unique name for the facility.
    T_min_app : float
        Minimum approach temperature observed during synthesis of heat exchanger network.

    Notes
    -----
    Original system Stream and HX objects are preserved. All Stream copies and new HX objects 
    can be found in a newly created flowsheet '<sys>_HXN' where <sys> is the system in which
    an instance of HeatExchangerNetwork is created.
    
    References
    ----------
    .. [1] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
        & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
        Heat Exchanger Networks (Chapter 9)
    
    Examples
    --------
    >>> from biosteam.units import ShortcutColumn, HXutility, Flash
    >>> from biosteam import Flowsheet
    >>> from biosteam.units.facilities import HeatExchangerNetwork
    >>> from biosteam import Stream, settings
    >>> from biosteam import main_flowsheet as f
    >>> flowsheet = Flowsheet('trial')
    >>> f.set_flowsheet(flowsheet)
    >>> settings.set_thermo(['Water', 'Methanol', 'Glycerol'])
    >>> feed1 = Stream('feed1', flow=(8000, 100, 25))
    >>> feed2 = Stream('feed2', flow=(10000, 1000, 10))
    >>> D1 = ShortcutColumn('D1', ins=feed1,
    ...                     outs=('distillate', 'bottoms_product'),
    ...                     LHK=('Methanol', 'Water'),
    ...                     y_top=0.99, x_bot=0.01, k=2,
    ...                     is_divided=True)
    >>> D1_H1 = HXutility('D1_H1', ins = D1.outs[1], T = 300)
    >>> D1_H2 = HXutility('D1_H2', ins = D1.outs[0], T = 300)
    >>> F1 = Flash('F1', ins=feed2,
    ...                     outs=('vapor', 'liquid'), V = 0.9, P = 101325)
    >>> HXN = HeatExchangerNetwork('trial_HXN', T_min_app = 5.)
    >>> trial_sys = f.create_system(
    ... 'trial_sys', feeds=[i for i in f.stream
    ...                  if i.sink and not i.source])
    >>> trial_sys.simulate()
    >>> HXN.simulate()
    >>> # See all results
    >>> round(HXN.actual_heat_util_load/HXN.original_heat_util_load, 2)
    0.88
    >>> HXN.get_stream_life_cycles()
    [<StreamLifeCycle: Stream_0, cold
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_4_hs>, H_in = 5.91e+06 kJ, H_out = 6.64e+06 kJ>
     		<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 6.64e+06 kJ, H_out = 4.26e+07 kJ>
     		<LifeStage: <HXutility: Util_0_hs>, H_in = 4.26e+07 kJ, H_out = 7.11e+07 kJ>
     	]>,
     <StreamLifeCycle: Stream_1, cold
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 3.67e-06 kJ, H_out = 1.45e+04 kJ>
     		<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 1.45e+04 kJ, H_out = 6.21e+06 kJ>
     		<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 6.21e+06 kJ, H_out = 2.54e+07 kJ>
     		<LifeStage: <HXutility: Util_1_hs>, H_in = 2.54e+07 kJ, H_out = 4.6e+08 kJ>
     	]>,
     <StreamLifeCycle: Stream_2, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 4.52e+07 kJ, H_out = 9.28e+06 kJ>
     		<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 9.28e+06 kJ, H_out = 3.08e+06 kJ>
     		<LifeStage: <HXutility: Util_2_hs>, H_in = 3.08e+06 kJ, H_out = 1.14e+06 kJ>
     	]>,
     <StreamLifeCycle: Stream_3, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 2.18e+07 kJ, H_out = 2.6e+06 kJ>
     	]>,
     <StreamLifeCycle: Stream_4, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_4_hs>, H_in = 7.49e+05 kJ, H_out = 2.24e+04 kJ>
     		<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 2.24e+04 kJ, H_out = 7.91e+03 kJ>
     		<LifeStage: <HXutility: Util_4_hs>, H_in = 7.91e+03 kJ, H_out = 2.91e+03 kJ>
     	]>]
        
    """

    network_priority = -1
    _N_ins = 0
    _N_outs = 0
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', T_min_app=5.):
        Facility.__init__(self, ID, None, None)
        self.T_min_app = T_min_app
        
    def _run(self):
        pass
    
    def _cost(self):
        sys = self.system
        sysname = sys.ID
        hx_utils = bst.process_tools.heat_exchanger_utilities_from_units(sys.units)
        hx_utils = [i for i in hx_utils if i.duty]
        hx_utils.sort(key = lambda x: x.duty)
        matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, actual_heat_util_load,\
        actual_cool_util_load, HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
        T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams, stream_HXs_dict,\
        hot_indices, cold_indices, original_heat_util_load, original_cool_util_load =\
        synthesize_network(hx_utils, ID_original=sysname, T_min_app=self.T_min_app)
        original_purchase_costs= [hx.purchase_cost for hx in hxs]
        original_installed_costs = [hx.installed_cost for hx in hxs]
        new_purchase_costs_HXp = []
        new_purchase_costs_HXu = []
        new_installed_costs_HXp = []
        new_installed_costs_HXu = []
        new_utility_costs = []
        for hx in new_HX_utils:
            new_installed_costs_HXu.append(hx.installed_cost)
            new_purchase_costs_HXu.append(hx.purchase_cost)
            new_utility_costs.append(hx.utility_cost)
        new_HXs = HXs_hot_side + HXs_cold_side
        for new_HX in new_HXs:
            new_purchase_costs_HXp.append(new_HX.purchase_cost)
            new_installed_costs_HXp.append(new_HX.installed_cost)
        
        self.cold_indices = cold_indices
        self.new_HXs = new_HXs
        self.new_HX_utils = new_HX_utils
        self.streams = streams
        
        stream_life_cycles = self.stream_life_cycles = self.get_stream_life_cycles()
        for life_cycle in stream_life_cycles:
            stages = life_cycle.life_cycle
            stream = life_cycle.index
            if stream not in cold_indices:
                for i in range(len(stages) - 1):
                    stage = stages[i]
                    hxn_unit = stage.unit
                    if hxn_unit.ID == 'Util_%s_hs'%stream:
                        new_HX_utils.remove(hxn_unit)
                        stream_in_at_stage = hxn_unit.ins[0]
                        
                        next_stage_unit = stages[i+1].unit
                        next_stage_unit_ID = next_stage_unit.ID
                        pointer = 1
                        if 'HX_%s_'%stream in next_stage_unit_ID\
                            or 'Util_%s_'%stream in next_stage_unit_ID:
                            pointer = 0
                        # try:
                        #     next_stage_unit.ins[pointer].vle(stream_in_at_stage.T, stream_in_at_stage.P)
                        #     next_stage_unit.simulate()
                        #     stages.remove(stage)
                        # except:
                        next_stage_unit.ins[pointer].T = stream_in_at_stage.T
                        # next_stage_unit.ins[pointer].H = stream_in_at_stage.H
                        next_stage_unit.simulate()
                        stages.remove(stage)
                        break
        self.purchase_costs['Heat exchangers'] = (sum(new_purchase_costs_HXp) + sum(new_purchase_costs_HXu)) \
            - (sum(original_purchase_costs))
        hu_sums1 = bst.HeatUtility.sum_by_agent(hx_utils_rearranged)
        hu_sums2 = bst.HeatUtility.sum_by_agent(sum([hx.heat_utilities for hx in new_HX_utils], ()))
        # to change sign on duty without switching heat/cool (i.e. negative costs):
        for hu in hu_sums1: hu.reverse()
        hus_final = tuple(bst.HeatUtility.sum_by_agent(hu_sums1 + hu_sums2))

        self._installed_cost = (sum(new_installed_costs_HXp) + sum(new_installed_costs_HXu)) \
            - (sum(original_installed_costs))
        self.heat_utilities = hus_final

        self.original_heat_utils = hx_utils_rearranged
        self.original_purchase_costs = original_purchase_costs
        self.original_utility_costs = hu_sums1
        self.new_purchase_costs_HXp = new_purchase_costs_HXp
        self.new_purchase_costs_HXu = new_purchase_costs_HXu
        self.new_utility_costs = hu_sums2
        self.stream_HXs_dict = stream_HXs_dict
        self.pinch_Ts = pinch_T_arr
        self.inlet_Ts = T_in_arr
        self.outlet_Ts = T_out_arr
        self.original_heat_util_load = original_heat_util_load
        self.original_cool_util_load = original_cool_util_load
        self.actual_heat_util_load = actual_heat_util_load
        self.actual_cool_util_load = actual_cool_util_load
        
    @property
    def installed_cost(self):
        return self._installed_cost
    
    def _design(self): pass
    
    def get_stream_life_cycles(self):
        # if hasattr(self, 'stream_life_cycles'): 
        #     return self.stream_life_cycles
        cold_indices = self.cold_indices
        new_HXs = self.new_HXs
        new_HX_utils = self.new_HX_utils
        streams = self.streams
        indices = [i for i in range(len(streams))]
        SLCs = [StreamLifeCycle(index, index in cold_indices) for index in indices]
        for SLC in SLCs:
            SLC.get_life_cycle(new_HXs, new_HX_utils)
        stream_life_cycles = SLCs
        self.stream_life_cycles = stream_life_cycles
        return stream_life_cycles
        
    def get_original_hxs_associated_with_streams(self): # pragma: no cover
        original_units = self.system.units
        original_heat_utils = self.original_heat_utils
        original_hx_utils = [i.heat_exchanger for i in original_heat_utils]
        original_hxs = {}
        stream_index = 0
        for hx in original_hx_utils:
            if '.' in hx.ID: # Names like 'U.1', i.e. non-explicitly named unit (e.g. auxillary HX)
                for unit in original_units:
                    if isinstance(unit, bst.units.MultiEffectEvaporator):
                        for key, component in unit.components.items():
                            if isinstance(component, list):
                                for subcomponent in component:
                                    if subcomponent is hx:
                                         original_hxs[stream_index] = (unit, key)
                            elif component is hx:
                                original_hxs[stream_index] = (unit, key)
                    elif isinstance(unit, bst.units.BinaryDistillation)\
                        or isinstance(unit, bst.units.ShortcutColumn):
                        if unit.boiler is hx:
                            original_hxs[stream_index] = (unit, 'boiler')
                        elif unit.condenser is hx:
                            original_hxs[stream_index] = (unit, 'condenser')
                    elif hasattr(unit, 'heat_exchanger'):
                        if unit.heat_exchanger is hx:
                            original_hxs[stream_index] = (unit, 'heat exchanger')
            else: # Explicitly named unit
                original_hxs[stream_index] = (hx, '')
            stream_index += 1
        self.original_hxs = original_hxs
        return original_hxs
    
    def save_stream_life_cycles_as_csv(self): # pragma: no cover
        if not hasattr(self, 'stream_life_cycles'):
            self.stream_life_cycles = self.get_stream_life_cycles()
        stream_life_cycles = self.stream_life_cycles
        if not hasattr(self, 'original_hxs'):
            self.original_hxs = self.get_original_hxs_associated_with_streams()
        original_hxs = self.original_hxs
        import csv
        from datetime import datetime
        dateTimeObj = datetime.now()
        filename = 'HXN-%s_%s.%s.%s.%s.%s.csv'%(self.system.ID, dateTimeObj.year,
                                                dateTimeObj.month, dateTimeObj.day,
                                                dateTimeObj.hour, dateTimeObj.minute)
        csvWriter = csv.writer(open(filename, 'w'), delimiter=',')
        csvWriter.writerow(['Stream', 'Type', 'Original unit', 'HXN unit', 'H_in (kJ)',
                            'H_out (kJ)', 'T_in (C)', 'T_out (C)'])
        stream, streamtype, original_unit, hxn_unit, H_in, H_out, T_in, T_out =\
            0, 0, 0, 0, 0, 0, 0, 0
            
        inlet_Ts = self.inlet_Ts
        outlet_Ts = self.outlet_Ts
        for life_cycle in stream_life_cycles:
            stream = life_cycle.index
            streamtype = 'Cold' if life_cycle.cold else 'Hot'
            stage_no = 0
            stages = life_cycle.life_cycle
            len_stages = len(stages)
            for stage in stages:
                original_unit = original_hxs[stream][0].ID
                if original_hxs[stream][1] is not '':
                     original_unit+= ' - ' + original_hxs[stream][1]
                
                hxn_unit = stage.unit
                hxn_unit_ID = hxn_unit.ID
                H_in = stage.H_in
                H_out = stage.H_out
                T_in, T_out = None, None
                if stage_no == 0:
                    T_in = inlet_Ts[stream] - 273.15
                if stage_no == len_stages - 1:
                    T_out = outlet_Ts[stream] - 273.15
                    
                row = [stream, streamtype, original_unit, hxn_unit_ID,
                       H_in, H_out, T_in, T_out]
                csvWriter.writerow(row)
                stage_no += 1