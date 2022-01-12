# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:44:24 2020

@author: sarangbhagwat
"""
import numpy as np
import copy
import biosteam as bst

__all__ = ('StreamLifeCycle', 'synthesize_network')

class LifeStage:
        
    def __init__(self, unit, index):
        self.unit = unit
        self.index = index
    
    @property
    def s_in(self): return self.unit.ins[self.index]
    
    @property
    def s_out(self): return self.unit.outs[self.index]
    
    @property
    def H_in(self): return self.s_in.H
    
    @property
    def H_out(self): return self.s_out.H
    
    def _info(self, N_tabs=1):
        tabs = N_tabs*'\t'
        return (f"{type(self).__name__}: {self.unit.ID}\n"
                + tabs + f"H_in = {self.H_in:.3g} kJ\n"
                + tabs + f"H_out = {self.H_out:.3g} kJ")
                
    def __repr__(self):
        return (f"<{type(self).__name__}: {repr(self.unit)}, H_in = {round(self.H_in, 4):.3g} kJ, H_out = {round(self.H_out, 4):.3g} kJ>")
        
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
        new_HXs_relevant, new_HX_utils_relevant =\
            self.get_relevant_units(index, new_HXs, new_HX_utils)
        life_cycle = (
            [LifeStage(unit, 0) for unit in new_HXs_relevant if name + '_' in unit.ins[0].ID]
            + [LifeStage(unit, 1) for unit in new_HXs_relevant if name + '_' in unit.ins[1].ID]
            + [LifeStage(unit, 0) for unit in new_HX_utils_relevant if name + '_' in unit.ins[0].ID]
        )
        life_cycle.sort(key = lambda pt: pt.H_in, reverse = not cold)
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


class Working_Life_Cycle:
    
    def __init__(self, index, cold):
        self.index = index
        self.name = 's_%s'%index
        self.cold = cold
        self.life_cycle = life_cycle = {}
        life_cycle['cold_side'] = []
        life_cycle['hot_side'] = []
        
    def add_stage(self, s_in, s_out, side):
        self.life_cycle[side].append(LifeStage(s_in, s_out))
    
    def sort_stages(self):
        life_cycle = self.life_cycle
        reverse = not self.cold
        life_cycle['cold_side'].sort(key = lambda stage: stage.H, reverse = reverse)
        life_cycle['cold_side'].sort(key = lambda stage: stage.H, reverse = reverse)
        hs_life_cycle = life_cycle['hot_side']
    
    def get_sorted_life_cycle(self):
        self.sort_stages()
        return self.life_cycle
    

def temperature_interval_pinch_analysis(hus, T_min_app = 10, force_ideal_thermo=False):
    hx_utils = hus
    hus_heating = [hu for hu in hx_utils if hu.duty > 0]
    hus_cooling = [hu for hu in hx_utils if hu.duty < 0]
    hx_utils_rearranged = hus_heating + hus_cooling
    hxs = [hu.heat_exchanger for hu in hx_utils_rearranged]
    if force_ideal_thermo:
        streams_inlet = [hx.ins[0] for hx in hxs]
        streams_quenched = [i.outs[0] for i in hxs]
        streams_inlet = [i.copy(thermo=i.thermo.ideal()) for i in streams_inlet]
        streams_quenched = [i.copy(thermo=i.thermo.ideal()) for i in streams_quenched]
    else:
        streams_inlet = [hx.ins[0].copy() for hx in hxs]
        streams_quenched = [i.outs[0].copy() for i in hxs]
    for i in streams_quenched: i.vle(H=i.H, P=i.P)
    for i in range(len(streams_inlet)):
        stream = streams_inlet[i]
        ID = 'Util_%s'%i
        stream.ID = 's_%s__%s'%(i,ID)
    N_heating = len(hus_heating)
    is_cold_stream_index = lambda x: x < N_heating
    T_in_arr = np.array([stream.T for stream in streams_inlet])
    T_out_arr = np.array([i.T for i in streams_quenched])
    adj_T_in_arr = T_in_arr.copy()
    adj_T_in_arr[:N_heating] -= T_min_app
    adj_T_out_arr = T_out_arr.copy()
    adj_T_out_arr[:N_heating] -= T_min_app
    T_changes_tuples = list(zip(adj_T_in_arr, adj_T_out_arr))
    all_Ts_descending = [*adj_T_in_arr, *adj_T_out_arr]
    all_Ts_descending.sort(reverse=True)
    stream_indices_for_T_intervals =\
        {(all_Ts_descending[i], all_Ts_descending[i+1]):[]\
            for i in range(len(all_Ts_descending)-1)}
    H_for_T_intervals = dict.fromkeys(stream_indices_for_T_intervals, 0)
    cold_indices = list(range(N_heating))
    hot_indices = list(range(N_heating, len(hxs)))
    indices = cold_indices + hot_indices
    for i in range(len(all_Ts_descending)-1):
        T_start = all_Ts_descending[i]
        T_end = all_Ts_descending[i+1]
        for stream_index in range(len(T_changes_tuples)):
            T1, T2 = T_changes_tuples[stream_index]
            if (T1 >= T_start and T2 <= T_end) or (T2 >= T_start and T1 <= T_end):
                multiplier = -1 if is_cold_stream_index(stream_index) else 1
                stream = streams_inlet[stream_index].copy()
                if stream.T != T_start: stream.vle(T = T_start, P = stream.P)
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
    # assert hot_util_load>= 0, 'Hot utility load is negative'
    # print(hot_util_load)
    # the lower temperature of the temperature interval for which the res_H is minimum
    pinch_cold_stream_T = all_Ts_descending[res_H_vector.index(-hot_util_load)+1]
    pinch_hot_stream_T = pinch_cold_stream_T + T_min_app
    cold_util_load = res_H_vector[len(res_H_vector)-1] + hot_util_load
    assert cold_util_load>=0, 'Cold utility load is negative'
    pinch_T_arr = []
    for i in cold_indices:
        if T_in_arr[i] > pinch_cold_stream_T:
            pinch_T_arr.append(T_in_arr[i])
        elif T_out_arr[i] < pinch_cold_stream_T:
            pinch_T_arr.append(T_out_arr[i])
        else:
            pinch_T_arr.append(pinch_cold_stream_T)
    for i in hot_indices:
        if T_in_arr[i] < pinch_hot_stream_T:
            pinch_T_arr.append(T_in_arr[i])
        elif T_out_arr[i] > pinch_hot_stream_T:
            pinch_T_arr.append(T_out_arr[i])
        else:
            pinch_T_arr.append(pinch_hot_stream_T)
    pinch_T_arr = np.array(pinch_T_arr)
    return pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr,\
           hxs, hot_indices, cold_indices, indices, streams_inlet, hx_utils_rearranged, \
           streams_quenched
            
        
def load_duties(streams, streams_quenched, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side):
    for index in indices:
        stream = streams[index].copy()
        H_in = stream.H
        stream.vle(T = pinch_T_arr[index], P = stream.P)
        H_pinch = stream.H
        H_out = streams_quenched[index].H
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


def synthesize_network(hus, T_min_app=5., Qmin=1e-3, force_ideal_thermo=False):  
    pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr,\
        hxs, hot_indices, cold_indices, indices, streams_inlet, hx_utils_rearranged, \
        streams_quenched = temperature_interval_pinch_analysis(hus, T_min_app, force_ideal_thermo)        
    H_out_arr = [i.H for i in streams_quenched]
    duties = np.array([abs(hx.Q)  for hx in hxs])
    dTs = np.abs(T_in_arr - T_out_arr)
    dTs[dTs == 0.] = 1e-12
    C_flow_vector = duties/dTs
    Q_hot_side = {}
    Q_cold_side = {}
    stream_HXs_dict = {i:[] for i in indices}
    is_cold = lambda x: x in cold_indices
    load_duties(streams_inlet, streams_quenched, pinch_T_arr, T_out_arr, indices, is_cold, Q_hot_side, Q_cold_side)
    matches_hs = {i: [] for i in cold_indices}
    matches_cs = {i: [] for i in hot_indices}
    candidate_hot_streams = hot_indices.copy()
    candidate_cold_streams = cold_indices.copy()
    HXs_hot_side = []
    HXs_cold_side = []
    streams_transient_cold_side = [i.copy() for i in streams_inlet]
    streams_transient_hot_side = [i.copy() for i in streams_inlet]
    for i in hot_indices:
        s = streams_transient_cold_side[i]
        if s.T != pinch_T_arr[i]:
            s.vle(T=pinch_T_arr[i], P=s.P)
    for i in cold_indices:
        s = streams_transient_hot_side[i]
        if s.T != pinch_T_arr[i]:
            s.vle(T=pinch_T_arr[i], P=s.P)
    
    def get_stream_at_H_max(cold):
        s_cs = streams_transient_cold_side[cold]
        s_hs = streams_transient_hot_side[cold]
        return s_cs if s_cs.H > s_hs.H else s_hs
    
    def get_stream_at_H_min(hot):
        s_cs = streams_transient_cold_side[hot]
        s_hs = streams_transient_hot_side[hot]
        return s_cs if s_cs.H < s_hs.H else s_hs
    
    def get_T_transient_cold_side(index):
        return streams_transient_cold_side[index].T
    
    def get_T_transient_hot_side(index):
        return streams_transient_hot_side[index].T
    
    attempts = set()
    # ------------- Cold side design ------------- # 
    unavailables = set([i for i in hot_indices if T_out_arr[i] >= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_in_arr[i] >= pinch_T_arr[i]])
    for hot in hot_indices:
        stream_quenched = False
        potential_matches = []
        for cold in cold_indices:
            if (C_flow_vector[hot]>= C_flow_vector[cold] and
                    get_T_transient_cold_side(hot) > get_T_transient_cold_side(cold) + T_min_app and
                    (hot not in unavailables) and (cold not in unavailables) and
                    (cold not in matches_cs[hot]) and (cold in candidate_cold_streams)): 
                potential_matches.append(cold)
        potential_matches = sorted(
            potential_matches, 
            key = lambda pot_cold: min(C_flow_vector[hot], C_flow_vector[pot_cold])
                                    * (get_T_transient_cold_side(hot) 
                                       - get_T_transient_cold_side(pot_cold)
                                       - T_min_app),
            reverse = True
        )
        for cold in potential_matches:
            ID = 'HX_%s_%s_cs'%(hot, cold)
            if ID in attempts: continue
            attempts.add(ID)
            Q_hstr = Q_cold_side[hot][1]
            Q_cstr = Q_cold_side[cold][1]
            Q_res = Q_cstr - Q_hstr
            # if abs(get_T_transient_cold_side(cold) - pinch_T_arr[cold])<= 0.01:
            #     continue
            hot_stream = streams_transient_cold_side[hot].copy()
            cold_stream = streams_transient_cold_side[cold].copy()
            
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            cold_stream.ID = 's_%s__%s'%(cold,ID)
            hot_out = bst.Stream('%s__s_%s'%(ID,hot), thermo=hot_stream.thermo)
            cold_out = bst.Stream('%s__s_%s'%(ID,cold), thermo=cold_stream.thermo)
            H_lim = H_out_arr[hot]
            new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                     outs = (hot_out, cold_out), H_lim0 = H_lim, 
                     T_lim1 = pinch_T_arr[cold], dT = T_min_app,
                     thermo = hot_stream.thermo)
            new_HX._run()
            if abs(new_HX.Q )< Qmin: continue
            HXs_cold_side.append(new_HX)
            stream_HXs_dict[hot].append(new_HX)
            stream_HXs_dict[cold].append(new_HX)
            Q_cold_side[hot][1] -= new_HX.Q
            Q_cold_side[cold][1] -= new_HX.Q
            streams_transient_cold_side[hot] = new_HX.outs[0]
            streams_transient_cold_side[cold] = new_HX.outs[1]
            H_out = new_HX.outs[0].H
            assert H_out - new_HX.ins[0].H <= 0.
            stream_quenched = H_out < H_lim or np.allclose(H_out, H_lim)
            matches_cs[hot].append(cold)
            if stream_quenched:
                break
    
    # ------------- Hot side design ------------- #                            
    unavailables = set([i for i in hot_indices if T_in_arr[i] <= pinch_T_arr[i]])
    unavailables.update([i for i in cold_indices if T_out_arr[i] <= pinch_T_arr[i]])
    
    for cold in cold_indices:
        # streams_transient_hot_side[cold] = streams_transient_cold_side[cold].copy()
        # T_transient_hot_side[cold] = min(T_transient_hot_side[cold], get_T_transient_cold_side(cold))
        potential_matches = []
        for hot in candidate_hot_streams:
            if ((cold in matches_cs and hot in matches_cs[cold])
                or (cold in matches_hs and hot in matches_hs[cold])):
                break
            if (C_flow_vector[cold]>= C_flow_vector[hot] and
                    get_T_transient_hot_side(hot) > get_T_transient_hot_side(cold) + T_min_app and
                    (hot not in unavailables) and (cold not in unavailables) and
                    (hot not in matches_hs[cold]) and (hot in candidate_hot_streams)):
                potential_matches.append(hot)
                
        potential_matches = sorted(potential_matches, key = lambda x:
                                    (min(C_flow_vector[cold], C_flow_vector[x])
                                        * ( get_T_transient_hot_side(x)
                                          - get_T_transient_hot_side(cold) - T_min_app)),
                                    reverse = True)
        stream_quenched = False
        for hot in potential_matches:
            ID = 'HX_%s_%s_hs'%(cold, hot)
            if ID in attempts: continue
            attempts.add(ID)
            Q_hstr = Q_hot_side[hot][1]
            Q_cstr = Q_hot_side[cold][1]
            Q_res = Q_cstr - Q_hstr
            # if abs(get_T_transient_hot_side(hot) - pinch_T_arr[hot])< 1e-6:
            #     continue
            hot_stream = streams_transient_hot_side[hot].copy()
            cold_stream = streams_transient_hot_side[cold].copy()
            cold_stream.ID = 's_%s__%s'%(cold,ID)
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            hot_out = bst.Stream('%s__s_%s'%(ID,hot), thermo=hot_stream.thermo)
            cold_out = bst.Stream('%s__s_%s'%(ID,cold), thermo=cold_stream.thermo)
            H_lim = H_out_arr[cold]
            new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                     outs = (cold_out, hot_out), H_lim0 = H_lim,
                     T_lim1 = pinch_T_arr[hot], dT = T_min_app,
                     thermo = hot_stream.thermo)
            try: new_HX._run()
            except: continue
            if abs(new_HX.Q)< Qmin: continue
            HXs_hot_side.append(new_HX)
            stream_HXs_dict[hot].append(new_HX)
            stream_HXs_dict[cold].append(new_HX)
            Q_hot_side[hot][1] -= new_HX.Q
            Q_hot_side[cold][1] -= new_HX.Q
            streams_transient_hot_side[cold] = new_HX.outs[0]
            streams_transient_hot_side[hot] = new_HX.outs[1]
            H_out = new_HX.outs[0].H
            assert H_out - new_HX.ins[0].H >= 0.
            stream_quenched = H_out > H_lim or np.allclose(H_out, H_lim)
            matches_hs[cold].append(hot)
            if stream_quenched:
                break
    
    # Offset heating requirement on cold side
    for cold in cold_indices:
        if Q_cold_side[cold][0]=='heat' and Q_cold_side[cold][1]>0:
            for hot in hot_indices:
                ID = 'HX_%s_%s_cs'%(hot, cold)
                if ID in attempts: continue
                attempts.add(ID)
                T_cold_in = get_T_transient_cold_side(cold)
                T_hot_in = get_T_transient_cold_side(hot)
                # if ((cold in matches_cs and hot in matches_cs[cold])
                #     or (cold in matches_hs and hot in matches_hs[cold])):
                #     continue
                if (Q_cold_side[hot][0]=='cool' and Q_cold_side[hot][1]>0 and
                        T_hot_in - T_cold_in >= T_min_app):
                    # if abs(T_cold_in - pinch_T_arr[cold])<= 0.01:
                    #     continue
                    hot_stream = streams_transient_cold_side[hot].copy()
                    cold_stream = streams_transient_cold_side[cold].copy()
                    hot_stream.ID = 's_%s__%s'%(hot,ID)
                    cold_stream.ID = 's_%s__%s'%(cold,ID)
                    hot_out = bst.Stream('%s__s_%s'%(ID,hot), thermo=hot_stream.thermo)
                    cold_out = bst.Stream('%s__s_%s'%(ID,cold), thermo=cold_stream.thermo)
                    new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                             outs = (hot_out, cold_out), H_lim0 = H_out_arr[hot],
                             T_lim1 = T_out_arr[cold], dT = T_min_app,
                             thermo = hot_stream.thermo)
                    try: new_HX._run()
                    except: continue
                    if abs(new_HX.Q )< Qmin: continue
                    HXs_cold_side.append(new_HX)
                    stream_HXs_dict[hot].append(new_HX)
                    stream_HXs_dict[cold].append(new_HX)                    
                    Q_cold_side[hot][1] -= new_HX.Q
                    Q_cold_side[cold][1] -= new_HX.Q                 
                    streams_transient_cold_side[hot] = new_HX.outs[0]
                    streams_transient_cold_side[cold] = new_HX.outs[1]
                    matches_cs[hot].append(cold)
                 
    # Offset cooling requirement on hot side  
    for hot in hot_indices:
        stream_quenched = False
        if Q_hot_side[hot][0]=='cool' and Q_hot_side[hot][1]>0:
            for cold in cold_indices:
                ID = 'HX_%s_%s_hs'%(cold, hot)
                if ID in attempts: continue
                attempts.add(ID)
                # if ((hot in matches_cs and cold in matches_cs[hot])
                #     or (hot in matches_hs and cold in matches_hs[hot])):
                #     continue
                original_cold_stream = get_stream_at_H_max(cold)
                T_cold_in = original_cold_stream.T
                T_hot_in = get_T_transient_hot_side(hot)
                if (Q_hot_side[cold][0]=='heat' and Q_hot_side[cold][1]>0 and
                        T_hot_in - T_cold_in>= T_min_app): 
                    # if abs(T_hot_in - pinch_T_arr[hot])<= 0.01:
                    #     continue                        
                    cold_stream = original_cold_stream.copy()
                    hot_stream = streams_transient_hot_side[hot].copy()
                    cold_stream.ID = 's_%s__%s'%(cold,ID)
                    hot_stream.ID = 's_%s__%s'%(hot,ID)
                    hot_out = bst.Stream('%s__s_%s'%(ID,hot), thermo=hot_stream.thermo)
                    cold_out = bst.Stream('%s__s_%s'%(ID,cold), thermo=cold_stream.thermo)
                    H_lim = H_out_arr[cold]
                    new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                             outs = (cold_out, hot_out), H_lim0 = H_lim, 
                             T_lim1 = T_out_arr[hot], dT = T_min_app,
                             thermo = hot_stream.thermo)
                    try: new_HX._run()
                    except: continue
                    if abs(new_HX.Q )< Qmin: continue
                    HXs_hot_side.append(new_HX)                        
                    stream_HXs_dict[hot].append(new_HX)
                    stream_HXs_dict[cold].append(new_HX)                        
                    Q_hot_side[hot][1] -= new_HX.Q
                    Q_hot_side[cold][1] -= new_HX.Q       
                    streams_transient_hot_side[cold] = new_HX.outs[0]
                    streams_transient_hot_side[hot] = new_HX.outs[1]
                    H_out = new_HX.outs[0].H
                    assert H_out - new_HX.ins[0].H >= 0.
                    stream_quenched = H_out > H_lim or np.allclose(H_out, H_lim)                   
                    matches_hs[cold].append(hot)                    
                    if stream_quenched:
                        break
    
    def get_hottest_stream_from_life_cycle(cold):
        s_hottest = get_stream_at_H_max(cold)
        H_max = s_hottest.H
        chemicals = s_hottest.chemicals
        for u in HXs_cold_side + HXs_hot_side:
            for s in u.outs:
                if '_%s_'%(cold) in s.ID and '__s_%s'%(cold) in s.ID and s.chemicals is chemicals and np.allclose(s.mol, s_hottest.mol):
                    if s.H > H_max:
                        H_max = s.H
                        s_hottest = s
        return s_hottest
    
    def get_coldest_stream_from_life_cycle(hot):
        s_coldest = get_stream_at_H_min(hot)
        H_min = s_coldest.H
        chemicals = s_coldest.chemicals
        for u in HXs_cold_side + HXs_hot_side:
            for s in u.outs:
                if '_%s_'%(hot) in s.ID and '__s_%s'%(hot) in s.ID and s.chemicals is chemicals and np.allclose(s.mol, s_coldest.mol):
                    if s.H < H_min:
                        H_min = s.H
                        s_coldest = s
        return s_coldest
    
    
    # Add final utility HXs
    new_HX_utils = []    
    for hot in hot_indices:
        hot_stream = get_coldest_stream_from_life_cycle(hot).copy()
        ID = 'Util_%s_cs'%(hot)
        hot_stream.ID = 's_%s__%s'%(hot,ID)
        outsID = '%s__s_%s'%(ID,hot)
        new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream, outs = outsID,
                                          H = H_out_arr[hot], rigorous = True,
                                          thermo = hot_stream.thermo)
        new_HX_util._run()
        s_out = new_HX_util-0
        np.testing.assert_allclose(s_out.H, H_out_arr[hot], rtol=5e-3, atol=1.)
        atol_T = 2. if 's' in hxs[hot].outs[0].phases else 0.001
        np.testing.assert_allclose(s_out.T, T_out_arr[hot], rtol=5e-3, atol=atol_T)
        new_HX_utils.append(new_HX_util)
        stream_HXs_dict[hot].append(new_HX_util)
            
    for cold in cold_indices:
        cold_stream = get_hottest_stream_from_life_cycle(cold).copy()
        ID = 'Util_%s_hs'%(cold)
        cold_stream.ID = 's_%s__%s'%(cold,ID)
        outsID = '%s__s_%s'%(ID,cold)
        new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, outs = outsID,
                                          H = H_out_arr[cold], rigorous = True,
                                          thermo = cold_stream.thermo)
        new_HX_util._run()
        s_out = new_HX_util-0
        np.testing.assert_allclose(s_out.H, H_out_arr[cold], rtol=1e-2, atol=1.)
        atol_T = 2. if 's' in hxs[cold].outs[0].phases else 0.001
        np.testing.assert_allclose(s_out.T, T_out_arr[cold], rtol=5e-2, atol=atol_T)
        new_HX_utils.append(new_HX_util)
        stream_HXs_dict[cold].append(new_HX_util)
    
    return HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
           T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams_inlet, stream_HXs_dict,\
           hot_indices, cold_indices

