# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Created on Sat Aug 22 21:58:19 2020
@author: sarangbhagwat
"""
from biosteam import HeatUtility, Facility
import biosteam as bst
from hxn.hxn_synthesis import synthesize_network, StreamLifeCycle

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
    >>> from hxn._heat_exchanger_network import HeatExchangerNetwork
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
    >>> round(HXN.act_heat_util_load/HXN.orig_heat_util_load, 3)
    0.884
    >>> HXN.get_stream_life_cycles()
    [<StreamLifeCycle: Stream_0, cold
     	<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 6.52e+06 kJ, H_out = 4.26e+07 kJ>
     	<LifeStage: <HXutility: Util_0_hs>, H_in = 4.26e+07 kJ, H_out = 7.35e+07 kJ>>,
     <StreamLifeCycle: Stream_1, cold
     	<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 3.6e-06 kJ, H_out = 6.08e+06 kJ>
     	<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 6.08e+06 kJ, H_out = 6.85e+06 kJ>
     	<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 6.85e+06 kJ, H_out = 2.78e+07 kJ>
     	<LifeStage: <HXutility: Util_1_hs>, H_in = 2.78e+07 kJ, H_out = 4.84e+08 kJ>>,
     <StreamLifeCycle: Stream_2, hot
     	<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 4.52e+07 kJ, H_out = 9.15e+06 kJ>
     	<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 9.15e+06 kJ, H_out = 3.07e+06 kJ>
     	<LifeStage: <HXutility: Util_2_hs>, H_in = 3.07e+06 kJ, H_out = 1.14e+06 kJ>>,
     <StreamLifeCycle: Stream_3, hot
     	<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 2.37e+07 kJ, H_out = 2.71e+06 kJ>>,
     <StreamLifeCycle: Stream_4, hot
     	<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 7.84e+05 kJ, H_out = 1.96e+04 kJ>
     	<LifeStage: <HXutility: Util_4_hs>, H_in = 1.96e+04 kJ, H_out = 2.91e+03 kJ>>]
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
        hx_utils.sort(key = lambda x: x.duty)
        matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, act_heat_util_load,\
        act_cool_util_load, HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
        T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams, stream_HXs_dict,\
        hot_indices, cold_indices, orig_heat_util_load, orig_cool_util_load =\
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
        self.purchase_costs['Heat exchangers'] = (sum(new_purchase_costs_HXp) + sum(new_purchase_costs_HXu)) \
            - (sum(original_purchase_costs))
        hu_sums1 = HeatUtility.sum_by_agent(hx_utils_rearranged)
        hu_sums2 = HeatUtility.sum_by_agent(sum([hx.heat_utilities for hx in new_HX_utils], ()))
        # to change sign on duty without switching heat/cool (i.e. negative costs):
        for hu in hu_sums1: hu.reverse()
        hus_final = tuple(HeatUtility.sum_by_agent(hu_sums1 + hu_sums2))
        self._installed_cost = (sum(new_installed_costs_HXp) + sum(new_installed_costs_HXu)) \
            - (sum(original_installed_costs))
        self.heat_utilities = hus_final
        self.cold_indices = cold_indices
        self.new_HXs = new_HXs
        self.new_HX_utils = new_HX_utils
        self.orig_heat_utils = hx_utils_rearranged
        self.original_purchase_costs = original_purchase_costs
        self.original_utility_costs = hu_sums1
        self.new_purchase_costs_HXp = new_purchase_costs_HXp
        self.new_purchase_costs_HXu = new_purchase_costs_HXu
        self.new_utility_costs = hu_sums2
        self.stream_HXs_dict = stream_HXs_dict
        self.streams = streams
        self.orig_heat_util_load = orig_heat_util_load
        self.orig_cool_util_load = orig_cool_util_load
        self.act_heat_util_load = act_heat_util_load
        self.act_cool_util_load = act_cool_util_load
        
    @property
    def installed_cost(self):
        return self._installed_cost
    
    def _design(self): pass
    
    def get_stream_life_cycles(self):
        cold_indices = self.cold_indices
        new_HXs = self.new_HXs
        new_HX_utils = self.new_HX_utils
        indices = [i for i in range(len(self.streams))]
        SLCs = [StreamLifeCycle(index, index in cold_indices) for index in indices]
        for SLC in SLCs:
            SLC.get_life_cycle(new_HXs, new_HX_utils)
        stream_life_cycles = SLCs
        self.stream_life_cycles = stream_life_cycles
        return stream_life_cycles