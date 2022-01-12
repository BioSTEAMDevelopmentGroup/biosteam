# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
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
import numpy as np
from .hxn_synthesis import synthesize_network, StreamLifeCycle
from warnings import warn

__all__ = ('HeatExchangerNetwork',)


class HeatExchangerNetwork(Facility):
    """
    Create a HeatExchangerNetwork object that will perform a pinch analysis
    on the entire system's heating and cooling utility objects. The heat
    exchanger network reduces the heating and cooling utility requirements 
    of the system and may add additional capital cost.
    
    Parameters
    ----------
    ID : str
        Unique name for the facility.
    T_min_app : float
        Minimum approach temperature observed during synthesis of heat exchanger network.
    units : Iterable[Unit], optional
        All unit operations available to the heat exchanger network. Defaults
        to all unit operations in the system.
    
    Notes
    -----
    Original system stream and heat exchanger objects are preserved. All stream 
    copies and new HX objects can be found in a newly created flowsheet 
    '<sys>_HXN' where <sys> is the name of the system associated to the 
    HeatExchangerNetwork object.
    
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
    0.82
    >>> abs(HXN.energy_balance_percent_error) < 0.01
    True
    >>> HXN.stream_life_cycles
    [<StreamLifeCycle: Stream_0, cold
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 5.91e+06 kJ, H_out = 4.25e+07 kJ>
     		<LifeStage: <HXutility: Util_0_hs>, H_in = 4.25e+07 kJ, H_out = 7.11e+07 kJ>
     	]>,
     <StreamLifeCycle: Stream_1, cold
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 0 kJ, H_out = 3.35e+04 kJ>
     		<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 3.35e+04 kJ, H_out = 5.54e+06 kJ>
     		<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 5.54e+06 kJ, H_out = 2.47e+07 kJ>
     		<LifeStage: <HXutility: Util_1_hs>, H_in = 2.47e+07 kJ, H_out = 2.79e+08 kJ>
     	]>,
     <StreamLifeCycle: Stream_2, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_0_2_hs>, H_in = 4.52e+07 kJ, H_out = 8.6e+06 kJ>
     		<LifeStage: <HXprocess: HX_1_2_hs>, H_in = 8.6e+06 kJ, H_out = 3.1e+06 kJ>
     		<LifeStage: <HXutility: Util_2_cs>, H_in = 3.1e+06 kJ, H_out = 1.14e+06 kJ>
     	]>,
     <StreamLifeCycle: Stream_3, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_1_3_hs>, H_in = 2.18e+07 kJ, H_out = 2.6e+06 kJ>
     		<LifeStage: <HXutility: Util_3_cs>, H_in = 2.6e+06 kJ, H_out = 2.6e+06 kJ>
     	]>,
     <StreamLifeCycle: Stream_4, hot
     	life_cycle = [
     		<LifeStage: <HXprocess: HX_1_4_hs>, H_in = 7.49e+05 kJ, H_out = 7.15e+05 kJ>
     		<LifeStage: <HXutility: Util_4_cs>, H_in = 7.15e+05 kJ, H_out = 7.15e+05 kJ>
     	]>]
        
    """
    ticket_name = 'HXN'
    acceptable_energy_balance_error = 0.02
    raise_energy_balance_error = False
    network_priority = -1
    _N_ins = 0
    _N_outs = 0
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
              'Work': 'kW'}
    
    def __init__(self, ID='', T_min_app=5., units=None, ignored=None, Qmin=1e-3,
                 force_ideal_thermo=False, cache_network=False):
        Facility.__init__(self, ID, None, None)
        self.T_min_app = T_min_app
        self.units = units
        self.ignored = ignored
        self.Qmin = Qmin
        self.force_ideal_thermo = force_ideal_thermo
        self.cache_network = cache_network
        
    def _get_original_heat_utilties(self):
        sys = self.system
        if self.units:
            units = self.units
            if callable(units): units = units()
        else:
            units = sys.units
        ignored = self.ignored
        if ignored:
            if callable(ignored): ignored = ignored()
            ignored_hx_utils = sum([i.heat_utilities for i in ignored], ())
        else:
            ignored_hx_utils = ()
        hx_utils = bst.process_tools.heat_exchanger_utilities_from_units(units)
        return [i for i in hx_utils if i.duty and i not in ignored_hx_utils]
        
    def _run(self): pass
    def _design(self): pass
    
    def _cost(self):
        sys = self.system
        hx_utils = self._get_original_heat_utilties()
        original_flowsheet = sys.flowsheet
        HXN_ID = sys.ID + '_HXN'
        bst.main_flowsheet.set_flowsheet(HXN_ID)
        use_cached_network = False
        if self.cache_network and hasattr(self, 'original_heat_utils'):
            hxs_cache = self.original_heat_exchangers
            hxs = [hu.heat_exchanger for hu in hx_utils]
            hxs_dct = {(i.owner, i._ID): i for i in hxs}
            try:
                hxs = [hxs_dct[i.owner, i._ID] for i in hxs_cache]
            except: pass
            else: use_cached_network = len(hxs) == len(hx_utils)
        try:
            if use_cached_network:
                hx_utils_rearranged = [i.heat_utilities[0] for i in hxs]
                stream_life_cycles = self.stream_life_cycles
                new_HXs = self.new_HXs
                new_HX_utils = self.new_HX_utils
                for i, life_cycle in enumerate(stream_life_cycles):
                    hx = hxs[i]
                    s_util_in = hx.ins[0]
                    stage = life_cycle.life_cycle[0]
                    s_lc = stage.unit.ins[stage.index]
                    s_lc.copy_like(s_util_in)
                    s_util_out = hx.outs[0]
                    H = s_util_out.H
                    for lc in life_cycle.life_cycle:
                        if isinstance(lc.unit, bst.HXutility):
                            lc.unit.H = H
                        else:
                            setattr(lc.unit, f'H_lim{lc.index}', s_util_out.H)
                sys = self.HXN_sys
                for unit in sys.units:
                    for s_in, s_out in zip(unit.ins, unit.outs):
                        if isinstance(s_out, bst.MultiStream):
                            s_out.F_mol = s_in.F_mol
                            if not (s_out.mol == s_in.mol).all():
                                s_out.copy_flow(s_in)
                                s_out.vle(T=s_out.T, P=s_out.P)
                        else:
                            s_out.mol[:] = s_in.mol
                # for life_cycle in stream_life_cycles:
                #     heating = life_cycle.cold
                #     for lc in life_cycle.life_cycle:
                #         if isinstance(lc.unit, bst.HXutility):
                #             H = lc.unit.H
                #         else:
                #             H = getattr(lc.unit, f'H_lim{lc.index}')
                #         outlet = lc.unit.outs[lc.index]
                #         if heating:
                #             if outlet.H > H: outlet.H = H
                #         elif outlet.H < H:
                #             outlet.H = H
            else:
                hx_utils.sort(key = lambda x: x.duty)
                self.HXN_flowsheet = HXN_F = bst.main_flowsheet
                for i in HXN_F.registries: i.clear()
                HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
                T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams_inlet, stream_HXs_dict,\
                hot_indices, cold_indices = \
                synthesize_network(hx_utils, self.T_min_app, self.Qmin, self.force_ideal_thermo)
                new_HXs = HXs_hot_side + HXs_cold_side
                self.cold_indices = cold_indices
                self.original_heat_exchangers = hxs
                self.new_HXs = new_HXs
                self.new_HX_utils = new_HX_utils
                self.streams_inlet = streams_inlet
                stream_life_cycles = self._get_stream_life_cycles()
                self.stream_HXs_dict = stream_HXs_dict
                self.pinch_Ts = pinch_T_arr
                self.inlet_Ts = T_in_arr
                self.outlet_Ts = T_out_arr
                all_units = new_HXs + new_HX_utils
                IDs = set([i.ID for i in all_units])
                assert len(all_units) == len(IDs)
                for i, life_cycle in enumerate(stream_life_cycles):
                    stage = life_cycle.life_cycle[0]
                    s_util = hx_utils_rearranged[i].heat_exchanger.ins[0]
                    s_lc = stage.unit.ins[stage.index]
                    s_lc.copy_like(s_util)
                for life_cycle in stream_life_cycles:
                    s_out = None
                    for i in life_cycle.life_cycle:
                        unit = i.unit
                        if s_out: unit.ins[i.index] = s_out
                        s_out = unit.outs[i.index]
                self.HXN_sys = sys = bst.System.from_units(None, all_units)
                sys.converge_method = 'fixedpoint'
                for i in sys.subsystems: i.converge_method = 'fixedpoint'
            
            original_purchase_costs = [hx.purchase_cost for hx in hxs]
            original_installed_costs = [hx.installed_cost for hx in hxs]
            # # Handle special case for heat exchanger crossing the pinch
            # for hx in new_HXs:
            #     if all([isinstance(i.sink, bst.HXutility) for i in hx.outs]):
            #         hx.Tlim1 = None
            #         hx.Hlim1 = hx.outs[1].sink.H
            try: 
                sys._converge()
            except:
                warning = RuntimeWarning('heat exchanger network was not able to converge')
                for i in sys.units: i._run()
                warn(warning)
            for i in sys.units:
                i._summary()
            for i in range(len(stream_life_cycles)):
                s_util = hx_utils_rearranged[i].heat_exchanger.outs[0]
                lc = stream_life_cycles[i].life_cycle[-1]
                s_lc = lc.unit.outs[lc.index]
                IDs = tuple([i.ID for i in s_util.available_chemicals])
                if use_cached_network:
                    try:
                        np.testing.assert_allclose(s_util.imol[IDs], s_lc.imol[IDs])
                        np.testing.assert_allclose(s_util.P, s_lc.P, rtol=1e-3, atol=0.1)
                        try:
                            np.testing.assert_allclose(s_util.H, s_lc.H, rtol=1e-3, atol=1.)
                        except:
                            lc.unit.simulate()
                            np.testing.assert_allclose(s_util.H, s_lc.H, rtol=1e-3, atol=1.)
                    except:
                        msg = ("heat exchanger network cache algorithm failed, cached network ignored")
                        warn(msg, RuntimeWarning, stacklevel=2)
                        del self.original_heat_utils
                        self._cost()
                        return
                else:
                    np.testing.assert_allclose(s_util.imol[IDs], s_lc.imol[IDs])
                    np.testing.assert_allclose(s_util.P, s_lc.P, rtol=1e-3, atol=0.1)
                    np.testing.assert_allclose(s_util.H, s_lc.H, rtol=1e-3, atol=1.)
            new_purchase_costs_HXp = []
            new_purchase_costs_HXu = []
            new_installed_costs_HXp = []
            new_installed_costs_HXu = []
            new_utility_costs = []
            for hx in new_HX_utils:
                new_installed_costs_HXu.append(hx.installed_cost)
                new_purchase_costs_HXu.append(hx.purchase_cost)
                new_utility_costs.append(hx.utility_cost)
            for new_HX in new_HXs:
                new_purchase_costs_HXp.append(new_HX.purchase_cost)
                new_installed_costs_HXp.append(new_HX.installed_cost)
            hu_sums1 = bst.HeatUtility.sum_by_agent(hx_utils_rearranged)
            new_heat_utils = sum([hx.heat_utilities for hx in new_HX_utils], ())
            hu_sums2 = bst.HeatUtility.sum_by_agent(new_heat_utils)
            # to change sign on duty without switching heat/cool (i.e. negative costs):
            for hu in hu_sums1: hu.reverse()
            hus_final = tuple(bst.HeatUtility.sum_by_agent(hu_sums1 + hu_sums2))
            Q_bal = (
                (2.*sum([abs(i.Q) for i in new_HXs])
                 + sum([abs(i.duty * i.agent.heat_transfer_efficiency) for i in hu_sums2]))
                / sum([abs(i.duty * i.agent.heat_transfer_efficiency) for i in hu_sums1])
            )
            energy_balance_error = Q_bal - 1
            self.energy_balance_percent_error = 100 * energy_balance_error
            self.installed_costs['Heat exchangers'] = (
                    sum(new_installed_costs_HXp)
                    + sum(new_installed_costs_HXu)
                    - sum(original_installed_costs)
            )
            self.purchase_costs['Heat exchangers'] = self.baseline_purchase_costs['Heat exchangers'] = (
                sum(new_purchase_costs_HXp) 
                + sum(new_purchase_costs_HXu)
                - sum(original_purchase_costs)
            )
            self.heat_utilities = hus_final
            self.original_heat_utils = hx_utils_rearranged
            self.original_purchase_costs = original_purchase_costs
            self.original_utility_costs = hu_sums1
            self.new_purchase_costs_HXp = new_purchase_costs_HXp
            self.new_purchase_costs_HXu = new_purchase_costs_HXu
            self.new_utility_costs = hu_sums2
            new_hus = bst.process_tools.heat_exchanger_utilities_from_units(new_HX_utils)
            hus_heating = [hu for hu in hx_utils if hu.duty > 0]
            hus_cooling = [hu for hu in hx_utils if hu.duty < 0]
            self.original_heat_util_load = sum([hu.duty for hu in hus_heating])
            self.original_cool_util_load = sum([abs(hu.duty) for hu in hus_cooling])
            self.actual_heat_util_load = sum([hu.duty for hu in new_hus if hu.duty>0])
            self.actual_cool_util_load = sum([abs(hu.duty) for hu in new_hus if hu.duty<0])
            if abs(energy_balance_error) > self.acceptable_energy_balance_error:
                msg = ("heat exchanger network energy balance is off by "
                      f"{energy_balance_error:.2%} (an absolute error greater "
                      f"than {self.acceptable_energy_balance_error:.2%})")
                if self.raise_energy_balance_error:
                    raise RuntimeError(msg)
                else:
                    warn(msg, RuntimeWarning, stacklevel=2)
        finally:
            bst.main_flowsheet.set_flowsheet(original_flowsheet)
    
    def _energy_balance_error_contributions(self):
        original_ignored = ignored = self.ignored
        if ignored and callable(ignored): ignored = ignored()
        energy_balance_errors = {}
        for hu in self._get_original_heat_utilties():
            self.ignored = ignored + [hu.heat_exchanger]
            if hasattr(hu.heat_exchanger, 'owner'):
                ID = hu.heat_exchanger.owner.ID, hu.heat_exchanger.ID
            else:
                ID = hu.heat_exchanger.ID
            try:
                self.simulate()
            except: 
                energy_balance_errors[ID] = (hu, None)
            else:
                energy_balance_errors[ID] = (hu, self.energy_balance_percent_error)
        self.ignored = original_ignored
        return energy_balance_errors
            
    def _get_stream_life_cycles(self):
        cold_indices = self.cold_indices
        new_HXs = self.new_HXs
        new_HX_utils = self.new_HX_utils
        streams = self.streams_inlet
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
            self.stream_life_cycles = self._get_stream_life_cycles()
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
                if original_hxs[stream][1]:
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
