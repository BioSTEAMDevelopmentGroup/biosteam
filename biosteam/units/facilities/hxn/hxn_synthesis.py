# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:44:24 2020

@author: sarangbhagwat
"""
import numpy as np
import copy
import biosteam as bst
from warnings import warn

__all__ = ('StreamLifeCycle', 'synthesize_network')

class LifeStage:
        
        def __init__(self, s_in, s_out):
            self.s_in = s_in
            self.s_out = s_out
        
        @property
        def H_in(self): return self.s_in.H
        
        @property
        def H_out(self): return self.s_out.H
        
        @property
        def unit(self): return self.s_in.sink
        
        def _info(self, N_tabs=1):
            tabs = N_tabs*'\t'
            return (f"{type(self).__name__}: {self.unit.ID}\n"
                    + tabs + f"H_in = {self.H_in:.3g} kJ\n"
                    + tabs + f"H_out = {self.H_out:.3g} kJ")
                    
        def __repr__(self):
            return (f"<{type(self).__name__}: {repr(self.unit)}, H_in = {self.H_in:.3g} kJ, H_out = {self.H_out:.3g} kJ>")
            
        def show(self):
            print(self._info())
        _ipython_display_ = show
        
        
class StreamLifeCycle:
    
    def __init__(self, index, cold):
        self.index = index
        self.name = 's_%s'%index
        self.cold = cold
        self.life_cycle = None
        
    def get_relevant_units(self, index, new_HXs, new_HX_utils):
        new_HXs_relevant = [hx for hx in new_HXs if '_%s_'%index in hx.ID]
        new_HX_utils_relevant = [hx for hx in new_HX_utils if '_%s_'%index in hx.ID]
        return new_HXs_relevant, new_HX_utils_relevant
        
    def get_life_cycle(self, new_HXs, new_HX_utils):
        index = self.index
        name = self.name
        cold = self.cold
        us_index = '_%s_'%index
        new_HXs_relevant, new_HX_utils_relevant =\
            self.get_relevant_units(index, new_HXs, new_HX_utils)
        life_cycle = [LifeStage(unit.ins[0], unit.outs[0])\
                      for unit in new_HXs_relevant if name in unit.ins[0].ID]\
                      + [LifeStage(unit.ins[1], unit.outs[1])\
                      for unit in new_HXs_relevant if name in unit.ins[1].ID]\
                      + [LifeStage(unit.ins[0], unit.outs[0])\
                      for unit in new_HX_utils_relevant if name in unit.ins[0].ID]
        life_cycle.sort(key = lambda pt: pt.H_out, reverse = not cold)
        self.life_cycle = life_cycle
        return life_cycle
        
    def __repr__(self):
        life_cycle = self.life_cycle
        cold = self.cold
        if not self.life_cycle:
            return 'Not initialized; run StreamLifeCycle.get_life_cycle or\
                  HX_Network.get_stream_life_cycles first.' 
        else:
            index = self.index
            name = 'Stream_%s'%index
            strtype = 'cold' if cold else 'hot'
            rep = ''
            for LifeStage in life_cycle:
                line = '\t\t' + repr(LifeStage) + '\n'
                rep += line
            rep = '<StreamLifeCycle: ' + name + ', ' + strtype  + '\n\tlife_cycle = [\n' +  rep[:-1] + '\n\t]>'
            return rep
        
    def show(self):
        info = repr(self).replace('[', '').replace(']', '').replace('life_cycle =', 'life_cycle:')
        print(info[1:-1])
        
    _ipython_display_ = show
    
    
def temperature_interval_pinch_analysis(hus, T_min_app = 10, ID_original = None):
    hx_utils = [hu for hu in hus\
                if abs(hu.heat_exchanger.ins[0].T - hu.heat_exchanger.outs[0].T)>0.01]
    hus_heating = [hu for hu in hx_utils if hu.duty > 0]
    hus_cooling = [hu for hu in hx_utils if hu.duty < 0]
    hxs_heating = [hu.heat_exchanger for hu in hx_utils if hu.duty > 0]
    hxs_cooling = [hu.heat_exchanger for hu in hx_utils if hu.duty < 0]
    streams_heating = [hx.ins[0].copy() for hx in hxs_heating]
    streams_cooling = [hx.ins[0].copy() for hx in hxs_cooling]
    init_heat_util, init_cool_util, init_Q = 0, 0, 0
    hx_utils_rearranged = hus_heating + hus_cooling
    hxs = hxs_heating + hxs_cooling
    streams = streams_heating + streams_cooling
    for i in range(len(streams)):
        stream = streams[i]
        hx = hxs[i]
        stream.vle(H = stream.H, P = stream.P)
        ID = 'Util_%s'%i
        stream.ID = 's_%s__%s'%(i,ID)
        
    len_hxs_heating = len(hxs_heating)
    is_cold_stream_index = lambda x: x<len_hxs_heating
    T_in_arr = np.array([stream.T for stream in streams])
    T_out_arr = np.array([hx.outs[0].T for hx in hxs])
    T_hot_side_arr = np.array([stream.T for stream in streams_heating]\
                               + [hx.outs[0].T for hx in hxs_cooling])
    T_cold_side_arr = np.array([hx.outs[0].T for hx in hxs_heating]\
                                + [stream.T for stream in streams_cooling])
    adj_T_in_arr = T_in_arr.copy()
    adj_T_in_arr[:len(hxs_heating)] -= T_min_app
    adj_T_out_arr = T_out_arr.copy()
    adj_T_out_arr[:len(hxs_heating)] -= T_min_app
    T_changes_tuples = list(zip(adj_T_in_arr, adj_T_out_arr))
    all_Ts_descending = [*adj_T_in_arr, *adj_T_out_arr]
    all_Ts_descending.sort(reverse=True)
    stream_indices_for_T_intervals =\
        {(all_Ts_descending[i], all_Ts_descending[i+1]):[]\
            for i in range(len(all_Ts_descending)-1)}
    H_for_T_intervals = dict.fromkeys(stream_indices_for_T_intervals, 0)
    cold_indices = list(range(len_hxs_heating))
    hot_indices = list(range(len_hxs_heating, len(hxs)))
    
    for i in range(len(all_Ts_descending)-1):
        T_start = all_Ts_descending[i]
        T_end = all_Ts_descending[i+1]
        dT = T_start - T_end
        for stream_index in range(len(T_changes_tuples)):
            if ((T_changes_tuples[stream_index][0]>= T_start and
                    T_changes_tuples[stream_index][1]<= T_end) or
                    (T_changes_tuples[stream_index][1]>= T_start and
                     T_changes_tuples[stream_index][0]<= T_end)):
                stream_indices_for_T_intervals[(T_start, T_end)].append(stream_index)
                multiplier = 1
                if is_cold_stream_index(stream_index):
                    multiplier = -1
                stream = streams[stream_index].copy()
                stream.vle(T = T_start, P = stream.P)
                H1 = stream.H
                stream.vle(T = T_end, P = stream.P)
                H2 = stream.H
                H = multiplier*(H1 - H2)
                H_for_T_intervals[(T_start, T_end)] += H
                
    res_H_vector = []
    prev_res_H = 0
    for interval, H in H_for_T_intervals.items():
        res_H_vector.append(prev_res_H + H)
        prev_res_H = res_H_vector[len(res_H_vector)-1]
    hot_util_load = - min(res_H_vector)
    assert hot_util_load>= 0
    # the lower temperature of the temperature interval for which the res_H is minimum
    pinch_cold_stream_T = all_Ts_descending[res_H_vector.index(-hot_util_load)+1]
    pinch_hot_stream_T = pinch_cold_stream_T + T_min_app
    cold_util_load = res_H_vector[len(res_H_vector)-1] + hot_util_load
    assert cold_util_load>=0
    pinch_T_arr = []
    
    for i in range(len(T_in_arr)):
        if not is_cold_stream_index(i):
            if T_in_arr[i]<pinch_hot_stream_T:
                pinch_T_arr.append(T_in_arr[i])
            elif T_out_arr[i]>pinch_hot_stream_T:
                pinch_T_arr.append(T_out_arr[i])
            else:
                pinch_T_arr.append(pinch_hot_stream_T)
        else:
            if T_in_arr[i]>pinch_cold_stream_T:
                pinch_T_arr.append(T_in_arr[i])
            elif T_out_arr[i]<pinch_cold_stream_T:
                pinch_T_arr.append(T_out_arr[i])
            else:
                pinch_T_arr.append(pinch_cold_stream_T)
    pinch_T_arr = np.array(pinch_T_arr)
    
    return pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr,\
           T_hot_side_arr, T_cold_side_arr, hus_heating,\
           hus_cooling, hxs_heating, hxs_cooling, hxs, hot_indices,\
           cold_indices, streams, hx_utils_rearranged
            
        
def load_duties(streams, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side):
    for index in indices:
        stream = streams[index].copy()
        H_in = stream.H
        stream.T = pinch_T_arr[index]
        stream.vle(T = pinch_T_arr[index], P = stream.P)
        H_pinch = stream.H
        stream.T = T_out_arr[index]
        stream.vle(T = T_out_arr[index], P = stream.P)
        H_out = stream.H
        if not is_cold(index):
            dH1 = abs(H_pinch - H_in)
            dH2 = abs(H_out - H_pinch)
            if abs(dH1)<0.01: dH1 = 0
            if abs(dH2)<0.01: dH2 = 0
            Q_hot_side[index] = ['cool', dH1]
            Q_cold_side[index] = ['cool', dH2]
        else:
            dH1 = H_out - H_pinch
            dH2 = H_pinch - H_in
            if abs(dH1)<0.01: dH1 = 0
            if abs(dH2)<0.01: dH2 = 0
            Q_hot_side[index] = ['heat', dH1]
            Q_cold_side[index] = ['heat', dH2]
            
            
def get_T_transient(pinch_T_arr, indices, T_in_arr):
    T_transient = pinch_T_arr.copy()
    T_transient[indices] = T_in_arr[indices]
    return T_transient


def synthesize_network(hus, ID_original, T_min_app=5.):  
    bst.main_flowsheet.set_flowsheet('%s_HXN'%ID_original)
    pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr,\
        T_hot_side_arr, T_cold_side_arr, hus_heating, hus_cooling, hxs_heating,\
        hxs_cooling, hxs, hot_indices, cold_indices, streams, hx_utils_rearranged =\
        temperature_interval_pinch_analysis(hus, T_min_app = T_min_app, ID_original = None)        
    duties = np.array([abs(hx.Q) for hx in hxs])
    C_flow_vector = duties/np.abs(T_in_arr - T_out_arr)
    Q_hot_side = {}
    Q_cold_side = {}
    T_transient_hot_side = get_T_transient(pinch_T_arr, hot_indices, T_in_arr)
    T_transient_cold_side = get_T_transient(pinch_T_arr, cold_indices, T_in_arr)
    indices = hot_indices + cold_indices
    stream_HXs_dict = {i:[] for i in indices}
    is_cold = lambda x: x in cold_indices
    load_duties(streams, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side)
    
    matches_hs = {i: [] for i in cold_indices}
    matches_cs = {i: [] for i in hot_indices}
    candidate_hot_streams = list(hot_indices)
    candidate_cold_streams = list(cold_indices)
    HXs_hot_side = []
    HXs_cold_side = []
    
                            # ------------- Cold side design ------------- # 
    unavailables = []
    unavailables = set([i for i in hot_indices if T_out_arr[i] >= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_in_arr[i] >= pinch_T_arr[i]])
    for hot in hot_indices:
        stream_quenched = False
        original_hot_stream = streams[hot]
        # T_transient_cold_side[hot]  = min(T_transient_cold_side[hot], T_transient_hot_side[hot] )
        potential_matches = []
        for cold in cold_indices:
            if (C_flow_vector[hot]>= C_flow_vector[cold] and
                    T_transient_cold_side[hot] > T_transient_cold_side[cold] + T_min_app and
                    (hot not in unavailables) and (cold not in unavailables) and
                    (cold not in matches_cs[hot]) and (cold in candidate_cold_streams)): 
                potential_matches.append(cold)
        
        potential_matches = sorted(potential_matches, key = lambda pot_cold:
                      (min(C_flow_vector[hot], C_flow_vector[pot_cold])
                       * (T_transient_cold_side[hot] 
                          - T_transient_cold_side[cold] - T_min_app)),
                      reverse = True)
        for cold in potential_matches:
            original_cold_stream = streams[cold]
            try:
                Q_hstr = Q_cold_side[hot][1]
                Q_cstr = Q_cold_side[cold][1]
                Q_res = Q_cstr - Q_hstr
                hot_stream = original_hot_stream.copy()
                hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
                cold_stream = original_cold_stream.copy()
                cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
                if abs(T_transient_cold_side[cold] - pinch_T_arr[cold])<= 0.01:
                    continue
                ID = 'HX_%s_%s_cs'%(hot, cold)
                hot_stream.ID = 's_%s__%s'%(hot,ID)
                cold_stream.ID = 's_%s__%s'%(cold,ID)
                outsIDs = ('%s__s_%s'%(ID,hot), '%s__s_%s'%(ID,cold))
                new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                         outs = outsIDs, T_lim0 = T_out_arr[hot],
                         T_lim1 = pinch_T_arr[cold], dT = T_min_app)
                new_HX.simulate()
                HXs_cold_side.append(new_HX)
                stream_HXs_dict[hot].append(new_HX)
                stream_HXs_dict[cold].append(new_HX)
                Q_cold_side[hot][1] -= new_HX.Q
                Q_cold_side[cold][1] -= new_HX.Q
                T_transient_cold_side[hot] = new_HX.outs[0].T
                T_transient_cold_side[cold] = new_HX.outs[1].T
                stream_quenched = T_transient_cold_side[hot] <= T_out_arr[hot]
                matches_cs[hot].append(cold)
            except:
                pass
            if stream_quenched:
                break
                
                            # ------------- Hot side design ------------- #                            
    unavailables = set([i for i in hot_indices if T_in_arr[i] <= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_out_arr[i] <= pinch_T_arr[i]])
    for cold in cold_indices:
        stream_quenched = False
        original_cold_stream = streams[cold]
        T_transient_hot_side[cold] = min(T_transient_hot_side[cold], T_transient_cold_side[cold])
        for hot in candidate_hot_streams:
            if ((cold in matches_cs.keys() and hot in matches_cs[cold])
                or (cold in matches_hs.keys() and hot in matches_hs[cold])):
                pass
            potential_matches = []
            if (C_flow_vector[cold]>= C_flow_vector[hot] and
                    T_transient_hot_side[hot] > T_transient_hot_side[cold] + T_min_app and
                    (hot not in unavailables) and (cold not in unavailables) and
                    (hot not in matches_hs[cold]) and (hot in candidate_hot_streams)):
                potential_matches.append(hot)
                
        potential_matches = sorted(potential_matches, key = lambda pot_hot:
                                   (min(C_flow_vector[cold], C_flow_vector[pot_hot])
                                       * ( T_transient_hot_side[pot_hot] 
                                          - T_transient_hot_side[cold] - T_min_app)),
                                   reverse = True)
            
        for hot in potential_matches:
            original_hot_stream = streams[hot]
            try:
                Q_hstr = Q_hot_side[hot][1]
                Q_cstr = Q_hot_side[cold][1]
                Q_res = Q_cstr - Q_hstr
                cold_stream = original_cold_stream.copy()
                cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)
                hot_stream = original_hot_stream.copy()
                hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)
                if abs(T_transient_hot_side[hot] - pinch_T_arr[hot])<= 0.01:
                    continue
                ID = 'HX_%s_%s_hs'%(cold, hot)
                cold_stream.ID = 's_%s__%s'%(cold,ID)
                hot_stream.ID = 's_%s__%s'%(hot,ID)
                outsIDs = ('%s__s_%s'%(ID,cold), '%s__s_%s'%(ID,hot))
                new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                         outs = outsIDs, T_lim0 = T_out_arr[cold],
                         T_lim1 = pinch_T_arr[hot], dT = T_min_app)
                new_HX.simulate()
                HXs_hot_side.append(new_HX)
                stream_HXs_dict[hot].append(new_HX)
                stream_HXs_dict[cold].append(new_HX)
                Q_hot_side[hot][1] -= new_HX.Q
                Q_hot_side[cold][1] -= new_HX.Q
                T_transient_hot_side[cold] = new_HX.outs[0].T
                T_transient_hot_side[hot] = new_HX.outs[1].T
                stream_quenched = T_transient_hot_side[cold] >= T_out_arr[cold]
                matches_hs[cold].append(hot)
            except:
                pass
            if stream_quenched:
                break
                
    # Offset heating requirement on cold side
    for cold in cold_indices:
        original_cold_stream = streams[cold]
        if Q_cold_side[cold][0]=='heat' and Q_cold_side[cold][1]>0:
            for hot in hot_indices:
                if ((cold in matches_cs.keys() and hot in matches_cs[cold])
                    or (cold in matches_hs.keys() and hot in matches_hs[cold])):
                    break
                original_hot_stream = streams[hot]
                if (Q_cold_side[hot][0]=='cool' and Q_cold_side[hot][1]>0 and
                        T_transient_cold_side[hot] - T_transient_cold_side[cold] >= T_min_app):
                    try:
                        hot_stream = original_hot_stream.copy()
                        hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
                        cold_stream = original_cold_stream.copy()
                        cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
                        if abs(T_transient_cold_side[cold] - pinch_T_arr[cold])<= 0.01:
                            continue
                        ID = 'HX_%s_%s_cs'%(hot, cold)
                        hot_stream.ID = 's_%s__%s'%(hot,ID)
                        cold_stream.ID = 's_%s__%s'%(cold,ID)
                        outsIDs = ('%s__s_%s'%(ID,hot), '%s__s_%s'%(ID,cold))
                        new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                                 outs = outsIDs, T_lim0 = T_out_arr[hot],
                                 T_lim1 = pinch_T_arr[cold], dT = T_min_app)
                        new_HX.simulate()
                        HXs_cold_side.append(new_HX)
                        stream_HXs_dict[hot].append(new_HX)
                        stream_HXs_dict[cold].append(new_HX)                    
                        Q_cold_side[hot][1] -= new_HX.Q
                        Q_cold_side[cold][1] -= new_HX.Q                        
                        T_transient_cold_side[hot] = new_HX.outs[0].T
                        T_transient_cold_side[cold] = new_HX.outs[1].T                        
                        matches_cs[hot].append(cold)
                    except:
                        pass     
                    
    # Offset cooling requirement on hot side    
    for hot in hot_indices:
        original_hot_stream = streams[hot]
        if Q_hot_side[hot][0]=='cool' and Q_hot_side[hot][1]>0:
            for cold in cold_indices:
                if ((hot in matches_cs.keys() and cold in matches_cs[hot])
                    or (hot in matches_hs.keys() and cold in matches_hs[hot])):
                    break
                original_cold_stream = streams[cold]
                if (Q_hot_side[cold][0]=='heat' and Q_hot_side[cold][1]>0 and
                        T_transient_hot_side[hot] - T_transient_hot_side[cold] >= T_min_app):                    
                    try:
                        cold_stream = original_cold_stream.copy()
                        cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)                        
                        hot_stream = original_hot_stream.copy()
                        hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)                        
                        if abs(T_transient_hot_side[hot] - pinch_T_arr[hot])<= 0.01:
                            continue                        
                        ID = 'HX_%s_%s_hs'%(cold, hot)
                        cold_stream.ID = 's_%s__%s'%(cold,ID)
                        hot_stream.ID = 's_%s__%s'%(hot,ID)
                        outsIDs = ('%s__s_%s'%(ID,cold), '%s__s_%s'%(ID,hot))
                        new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                                 outs = outsIDs, T_lim0 = T_out_arr[cold], 
                                 T_lim1 = pinch_T_arr[hot], dT = T_min_app)
                        new_HX.simulate()
                        HXs_hot_side.append(new_HX)                        
                        stream_HXs_dict[hot].append(new_HX)
                        stream_HXs_dict[cold].append(new_HX)                        
                        Q_hot_side[hot][1] -= new_HX.Q
                        Q_hot_side[cold][1] -= new_HX.Q                        
                        T_transient_hot_side[cold] = new_HX.outs[0].T
                        T_transient_hot_side[hot] = new_HX.outs[1].T       
                        stream_quenched = T_transient_hot_side[cold] >= T_out_arr[cold]                    
                        matches_hs[cold].append(hot)                    
                    except:
                        pass     
                    
    # Add final utility HXs
    new_HX_utils = []    
    for hot in hot_indices:
        new_HX_util = None
        if T_transient_cold_side[hot] > T_out_arr[hot]:
            original_hot_stream = streams[hot]
            hot_stream = original_hot_stream.copy()
            hot_stream.vle(T = T_transient_cold_side[hot], P = hot_stream.P)
            ID = 'Util_%s_cs'%(hot)
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            outsID = '%s__s_%s'%(ID,hot)
            new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream,outs = outsID,
                                              T = T_out_arr[hot], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[hot].append(new_HX_util)
        if T_transient_hot_side[hot] > pinch_T_arr[hot] + 0.05:
            original_hot_stream = streams[hot]
            hot_stream = original_hot_stream.copy()
            hot_stream.vle(T = T_transient_hot_side[hot], P = hot_stream.P)
            ID = 'Util_%s_hs'%(hot)
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            outsID = '%s__s_%s'%(ID,hot)
            new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream, outs = outsID,
                                              T = pinch_T_arr[hot], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[hot].append(new_HX_util)
    for cold in cold_indices:
        if T_transient_hot_side[cold] < T_out_arr[cold]:
            original_cold_stream = streams[cold]
            cold_stream = original_cold_stream.copy()
            cold_stream.vle(T = T_transient_hot_side[cold], P = cold_stream.P)
            ID = 'Util_%s_hs'%(cold)
            cold_stream.ID = 's_%s__%s'%(cold,ID)
            outsID = '%s__s_%s'%(ID,cold)
            new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, outs = outsID,
                                              T = T_out_arr[cold], rigorous = True)
            new_HX_util.simulate()
            new_HX_utils.append(new_HX_util)
            stream_HXs_dict[cold].append(new_HX_util)
        # if T_transient_cold_side[cold] + 0.05 < pinch_T_arr[cold]:
        #     original_cold_stream = streams[cold]
        #     cold_stream = original_cold_stream.copy()
        #     cold_stream.vle(T = T_transient_cold_side[cold], P = cold_stream.P)
        #     ID = 'Util_%s_cs'%(cold)
        #     cold_stream.ID = 's_%s__%s'%(cold,ID)
        #     outsID = '%s__s_%s'%(ID,cold)
        #     new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, outs = outsID,
        #                                       T = pinch_T_arr[cold], rigorous = True)
        #     new_HX_util.simulate()
        #     new_HX_utils.append(new_HX_util)
        #     stream_HXs_dict[cold].append(new_HX_util)
            
    new_hus = bst.process_tools.heat_exchanger_utilities_from_units(new_HX_utils)
    act_cool_util_load = sum([abs(hu.duty) for hu in new_hus if hu.duty<0])
    act_heat_util_load = sum([hu.duty for hu in new_hus if hu.duty>0])
    orig_heat_util_load = sum([hu.duty for hu in hus_heating])
    orig_cool_util_load = sum([abs(hu.duty) for hu in hus_cooling])
    Q_prev_heating = sum([hx.Q for hx in hxs_heating])
    Q_prev_cooling = sum([abs(hx.Q) for hx in hxs_cooling])
    Q_prev = Q_prev_heating + Q_prev_cooling
    Q_HXp = sum([hx.Q for hx in HXs_hot_side]) + sum([hx.Q for hx in HXs_cold_side])
    Q_new_heating = sum([hx_util.Q for hx_util in new_HX_utils
                         if hx_util.ins[0].T < hx_util.outs[0].T])
    Q_new_cooling = sum([abs(hx_util.Q) for hx_util in new_HX_utils
                         if hx_util.ins[0].T > hx_util.outs[0].T])
    Q_new_utils = Q_new_heating + Q_new_cooling
    Q_new = 2*Q_HXp + Q_new_utils
    Q_bal = Q_new/Q_prev
    if abs(Q_bal - 1)>0.02:
        msg = '\n\n\n WARNING: Q balance of HXN off by %s p.c.,\ which is more than 2 p.c.\n\n\n'\
              %(100*(Q_bal - 1))
        warn(msg, UserWarning, stacklevel =2)
    
    return matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, act_heat_util_load,\
           act_cool_util_load, HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
           T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams, stream_HXs_dict,\
           hot_indices, cold_indices, orig_heat_util_load, orig_cool_util_load
