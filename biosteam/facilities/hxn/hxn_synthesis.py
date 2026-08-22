# -*- coding: utf-8 -*-
# HXN: The automated Heat Exchanger Network design package.
# Copyright (C) 2020-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/sarangbhagwat/hxn/blob/master/LICENSE.txt
# for license details.
"""
Created on Sat May  2 16:44:24 2020

@author: sarangbhagwat
"""
from collections import namedtuple
import numpy as np
import biosteam as bst
from warnings import warn

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
    
    def get_sorted_life_cycle(self):
        self.sort_stages()
        return self.life_cycle
    

ProblemTable = namedtuple(
    'ProblemTable',
    ['Ts', 'interval_H', 'point_H', 'residual',
     'hot_util_load', 'cold_util_load', 'pinch_T']
)

def _stream_H_at_boundaries(stream_in, H_in, H_out, T_lo, T_hi, Ts, shift):
    """
    Enthalpies [kJ/hr] of one monotone stream at the grid boundaries
    `Ts` (shifted scale, descending, all within [T_lo, T_hi]).

    Exact at the stream's own end points (H_in/H_out as given); in between,
    the inlet copy is flashed at the *real* temperature `T + shift` and the
    result is clipped to [min(H_in, H_out), max(H_in, H_out)] so that a
    non-equilibrium outlet (e.g. a column reboiler/condenser product) can
    never inflate an interval. A single copy is walked down the grid so
    each VLE is warm-started from the previous boundary.
    """
    H_lo, H_hi = sorted((H_in, H_out))
    H_top, H_bottom = (H_in, H_out) if H_in > H_out else (H_out, H_in)
    Hs = np.empty(Ts.size)
    stream = stream_in.copy()
    for k, T in enumerate(Ts):
        if T == T_hi:
            Hs[k] = H_top
        elif T == T_lo:
            Hs[k] = H_bottom
        else:
            T_real = T + shift
            try:
                stream.vle(T=T_real, P=stream.P)
                H = stream.H
            except Exception as error:
                warn(f"could not solve VLE for {stream!r} at {T_real:.2f} K "
                     f"({error!r}); interpolating enthalpy linearly in "
                     "temperature for the problem table", RuntimeWarning)
                H = H_lo + (H_hi - H_lo) * (T - T_lo) / (T_hi - T_lo)
            Hs[k] = min(max(H, H_lo), H_hi)
    return Hs

def problem_table(streams_inlet, streams_quenched, is_hot, T_min_app):
    """
    Energy-consistent problem table (temperature-interval heat cascade).

    Parameters
    ----------
    streams_inlet : list[Stream]
        Inlet stream of each utility heat exchanger.
    streams_quenched : list[Stream]
        Corresponding outlet streams, re-flashed at their enthalpy.
    is_hot : Sequence[bool]
        True where the stream is cooled.
    T_min_app : float
        Minimum approach temperature [K].

    Returns
    -------
    ProblemTable
        Grid temperatures `Ts` (shifted scale, descending), per-stream
        `interval_H` (N x n-1) and `point_H` (N x n) contributions (+ for
        hot, - for cold), the cascade `residual` (n), `hot_util_load`,
        `cold_util_load` and the shifted-scale `pinch_T`.

    Notes
    -----
    Hot streams are shifted down by `T_min_app`; cold streams are not. For
    monotone streams the contribution to interval (Ts[k], Ts[k+1]) is
    sign * (H(Ts[k]) - H(Ts[k+1])) with H evaluated at the real temperature
    and clipped to [H_in, H_out], so every stream's contributions telescope
    exactly to sign * |H_out - H_in|. Isothermal streams, and streams whose
    outlet temperature moves against their duty (a heated stream that exits
    colder than it entered, e.g. a reboiler outlet at VLE), are point loads
    at their outlet temperature. The cascade starting from zero hot utility
    is residual[k] = sum(point_H[:, :k+1]) + sum(interval_H[:, :k]); the
    minimum fixes the hot utility target, `residual[-1] + hot_util_load`
    the cold one, and its location the pinch. With the per-stream identity
    above, hot_util_load - cold_util_load equals the net heating demand.
    """
    N = len(streams_inlet)
    is_hot = np.asarray(is_hot, dtype=bool)
    sign = np.where(is_hot, 1., -1.)
    shift = np.where(is_hot, T_min_app, 0.)
    T_in = np.array([s.T for s in streams_inlet])
    T_out = np.array([s.T for s in streams_quenched])
    H_in = np.array([s.H for s in streams_inlet])
    H_out = np.array([s.H for s in streams_quenched])
    monotone = (sign * (T_in - T_out)) > 0.
    T_hi = np.where(monotone, np.maximum(T_in, T_out), T_out) - shift
    T_lo = np.where(monotone, np.minimum(T_in, T_out), T_out) - shift
    Ts = np.unique(np.concatenate([T_hi, T_lo]))[::-1]
    n = Ts.size
    interval_H = np.zeros((N, n - 1))
    point_H = np.zeros((N, n))
    for j in range(N):
        if monotone[j]:
            idx = np.flatnonzero((Ts <= T_hi[j]) & (Ts >= T_lo[j]))
            Hs = _stream_H_at_boundaries(streams_inlet[j], H_in[j], H_out[j],
                                         T_lo[j], T_hi[j], Ts[idx], shift[j])
            interval_H[j, idx[:-1]] = sign[j] * (Hs[:-1] - Hs[1:])
        else:
            k = np.searchsorted(-Ts, -T_hi[j])
            point_H[j, k] = sign[j] * abs(H_out[j] - H_in[j])
    residual = np.cumsum(
        point_H.sum(axis=0) + np.concatenate([[0.], interval_H.sum(axis=0)])
    )
    k_pinch = int(np.argmin(residual))
    scale = np.abs(H_out - H_in).sum()
    if -residual[k_pinch] <= 1e-9 * scale:  # threshold problem: no hot utility
        hot_util_load = 0.
        k_pinch = 0
    else:
        hot_util_load = -residual[k_pinch]
    cold_util_load = residual[-1] + hot_util_load
    return ProblemTable(Ts, interval_H, point_H, residual,
                        hot_util_load, cold_util_load, Ts[k_pinch])

def temperature_interval_pinch_analysis(hus,
                                        T_min_app=10,
                                        force_ideal_thermo=False,
                                        sort_hus_by_T=False):
    hx_utils = hus
    hus_heating = [hu for hu in hx_utils if hu.duty > 0]
    hus_cooling = [hu for hu in hx_utils if hu.duty < 0]
    if sort_hus_by_T:
        hus_heating.sort(key=lambda i: i.unit.ins[0].T, reverse=True)
        hus_cooling.sort(key=lambda i: i.unit.ins[0].T)
    hx_utils_rearranged = hus_heating + hus_cooling
    hxs = [hu.unit for hu in hx_utils_rearranged]
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
    cold_indices = list(range(N_heating))
    hot_indices = list(range(N_heating, len(hxs)))
    indices = cold_indices + hot_indices
    T_in_arr = np.array([stream.T for stream in streams_inlet])
    T_out_arr = np.array([i.T for i in streams_quenched])
    is_hot = np.zeros(len(hxs), dtype=bool)
    is_hot[hot_indices] = True
    table = problem_table(streams_inlet, streams_quenched, is_hot, T_min_app)
    hot_util_load = table.hot_util_load
    cold_util_load = table.cold_util_load
    pinch_cold_stream_T = table.pinch_T
    pinch_hot_stream_T = pinch_cold_stream_T + T_min_app
    # Per-stream pinch temperature: where the stream is split between the
    # hot-side and cold-side designs. Non-monotone streams (T_out against
    # the duty) are not split: their pinch is the inlet T, so load_duties
    # puts the whole duty on the hot side, as before this fix.
    pinch_T_arr = []
    for i in cold_indices:
        if T_in_arr[i] > pinch_cold_stream_T or T_in_arr[i] > T_out_arr[i]:
            pinch_T_arr.append(T_in_arr[i])
        elif T_out_arr[i] < pinch_cold_stream_T:
            pinch_T_arr.append(T_out_arr[i])
        else:
            pinch_T_arr.append(pinch_cold_stream_T)
    for i in hot_indices:
        if T_in_arr[i] < pinch_hot_stream_T or T_in_arr[i] < T_out_arr[i]:
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

def synthesize_network(hus, T_min_app=5., Qmin=1e-3, force_ideal_thermo=False,
                       avoid_recycle=False, sort_hus_by_T=False):  
    pinch_T_arr, hot_util_load, cold_util_load, T_in_arr, T_out_arr,\
        hxs, hot_indices, cold_indices, indices, streams_inlet, hx_utils_rearranged, \
        streams_quenched = temperature_interval_pinch_analysis(hus, T_min_app, force_ideal_thermo,
                                                               sort_hus_by_T)        
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
    HXs_hot_side = []
    HXs_cold_side = []
    streams_transient_cold_side = streams_inlet
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
    success = set()
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
                    (cold not in matches_cs[hot]) and (cold in cold_indices)): 
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
            match = (hot, cold)
            ID = 'HX_%s_%s_cs' % match
            if ID in attempts or (avoid_recycle and match in success): continue
            attempts.add(ID)
            hot_stream = streams_transient_cold_side[hot].copy()
            cold_stream = streams_transient_cold_side[cold].copy()
            
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            cold_stream.ID = 's_%s__%s'%(cold,ID)
            hot_out = hot_stream.copy('%s__s_%s'%(ID,hot))
            cold_out = cold_stream.copy('%s__s_%s'%(ID,cold))
            H_lim = H_out_arr[hot]
            new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                     outs = (hot_out, cold_out), H_lim0 = H_lim, 
                     T_lim1 = pinch_T_arr[cold], dT = T_min_app,
                     thermo = hot_stream.thermo)
            try: new_HX._run()
            except: continue
            if abs(new_HX.Q )< Qmin: continue
            success.add(match)
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
        potential_matches = []
        for hot in hot_indices:
            if ((cold in matches_cs and hot in matches_cs[cold])
                or (cold in matches_hs and hot in matches_hs[cold])):
                break
            if (C_flow_vector[cold]>= C_flow_vector[hot] and
                    get_T_transient_hot_side(hot) > get_T_transient_hot_side(cold) + T_min_app and
                    (hot not in unavailables) and (cold not in unavailables) and
                    (hot not in matches_hs[cold]) and (hot in hot_indices)):
                potential_matches.append(hot)
                
        potential_matches = sorted(potential_matches, key = lambda x:
                                    (min(C_flow_vector[cold], C_flow_vector[x])
                                        * ( get_T_transient_hot_side(x)
                                          - get_T_transient_hot_side(cold) - T_min_app)),
                                    reverse = True)
        stream_quenched = False
        for hot in potential_matches:
            match = (hot, cold)
            ID = 'HX_%s_%s_hs' % (cold, hot)
            if ID in attempts or (avoid_recycle and match in success): continue
            attempts.add(ID)
            hot_stream = streams_transient_hot_side[hot].copy()
            cold_stream = streams_transient_hot_side[cold].copy()
            cold_stream.ID = 's_%s__%s'%(cold,ID)
            hot_stream.ID = 's_%s__%s'%(hot,ID)
            hot_out = hot_stream.copy('%s__s_%s'%(ID,hot))
            cold_out = cold_stream.copy('%s__s_%s'%(ID,cold))
            H_lim = H_out_arr[cold]
            new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                     outs = (cold_out, hot_out), H_lim0 = H_lim,
                     T_lim1 = pinch_T_arr[hot], dT = T_min_app,
                     thermo = hot_stream.thermo)
            try: new_HX._run()
            except: continue
            if abs(new_HX.Q)< Qmin: continue
            success.add(match)
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
                match = (hot, cold)
                ID = 'HX_%s_%s_cs' % match
                if ID in attempts or (avoid_recycle and match in success): continue
                attempts.add(ID)
                T_cold_in = get_T_transient_cold_side(cold)
                T_hot_in = get_T_transient_cold_side(hot)
                if (Q_cold_side[hot][0]=='cool' and Q_cold_side[hot][1]>0 and
                        T_hot_in - T_cold_in >= T_min_app):
                    hot_stream = streams_transient_cold_side[hot].copy()
                    cold_stream = streams_transient_cold_side[cold].copy()
                    hot_stream.ID = 's_%s__%s'%(hot,ID)
                    cold_stream.ID = 's_%s__%s'%(cold,ID)
                    hot_out = hot_stream.copy('%s__s_%s'%(ID,hot))
                    cold_out = cold_stream.copy('%s__s_%s'%(ID,cold))
                    new_HX = bst.units.HXprocess(ID = ID, ins = (hot_stream, cold_stream),
                             outs = (hot_out, cold_out), H_lim0 = H_out_arr[hot],
                             T_lim1 = T_out_arr[cold], dT = T_min_app,
                             thermo = hot_stream.thermo)
                    try: new_HX._run()
                    except: continue
                    if abs(new_HX.Q )< Qmin: continue
                    success.add(match)
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
                match = (hot, cold)
                ID = 'HX_%s_%s_hs' % (cold, hot)
                if ID in attempts or (avoid_recycle and match in success): continue
                attempts.add(ID)
                original_cold_stream = get_stream_at_H_max(cold)
                T_cold_in = original_cold_stream.T
                T_hot_in = get_T_transient_hot_side(hot)
                if (Q_hot_side[cold][0]=='heat' and Q_hot_side[cold][1]>0 and
                        T_hot_in - T_cold_in>= T_min_app):    
                    cold_stream = original_cold_stream
                    hot_stream = streams_transient_hot_side[hot].copy()
                    cold_stream.ID = 's_%s__%s'%(cold,ID)
                    hot_stream.ID = 's_%s__%s'%(hot,ID)
                    hot_out = hot_stream.copy('%s__s_%s'%(ID,hot))
                    cold_out = cold_stream.copy('%s__s_%s'%(ID,cold))
                    H_lim = H_out_arr[cold]
                    new_HX = bst.units.HXprocess(ID = ID, ins = (cold_stream, hot_stream),
                             outs = (cold_out, hot_out), H_lim0 = H_lim, 
                             T_lim1 = T_out_arr[hot], dT = T_min_app,
                             thermo = hot_stream.thermo)
                    try: new_HX._run()
                    except: continue
                    if abs(new_HX.Q )< Qmin: continue
                    success.add(match)
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
    
    # Add final utility HXs
    new_HX_utils = []    
    for hot in hot_indices:
        hot_stream = get_stream_at_H_min(hot)
        ID = 'Util_%s_cs'%(hot)
        hot_stream.ID = 's_%s__%s'%(hot,ID)
        outlet = hot_stream.copy('%s__s_%s'%(ID,hot))
        new_HX_util = bst.units.HXutility(ID = ID, ins = hot_stream, outs = outlet,
                                          H = H_out_arr[hot], rigorous = True,
                                          thermo = hot_stream.thermo)
        new_HX_util._run()
        s_out = new_HX_util-0
        np.testing.assert_allclose(s_out.H, H_out_arr[hot], rtol=5e-3, atol=1.)
        atol_T = 5. if 's' in hxs[hot].outs[0].phases else 0.001
        np.testing.assert_allclose(s_out.T, T_out_arr[hot], rtol=5e-3, atol=atol_T)
        new_HX_utils.append(new_HX_util)
        stream_HXs_dict[hot].append(new_HX_util)
            
    for cold in cold_indices:
        cold_stream = get_stream_at_H_max(cold)
        ID = 'Util_%s_hs'%(cold)
        cold_stream.ID = 's_%s__%s'%(cold,ID)
        outlet = cold_stream.copy('%s__s_%s'%(ID,cold))
        new_HX_util = bst.units.HXutility(ID = ID, ins = cold_stream, outs = outlet,
                                          H = H_out_arr[cold], rigorous = True,
                                          thermo = cold_stream.thermo)
        new_HX_util._run()
        s_out = new_HX_util.outs[0]
        np.testing.assert_allclose(s_out.H, H_out_arr[cold], rtol=1e-2, atol=1.)
        atol_T = 5. if 's' in hxs[cold].outs[0].phases else 0.001
        np.testing.assert_allclose(s_out.T, T_out_arr[cold], rtol=5e-2, atol=atol_T)
        new_HX_utils.append(new_HX_util)
        stream_HXs_dict[cold].append(new_HX_util)
    
    return HXs_hot_side, HXs_cold_side, new_HX_utils, hxs, T_in_arr,\
           T_out_arr, pinch_T_arr, C_flow_vector, hx_utils_rearranged, streams_inlet, stream_HXs_dict,\
           hot_indices, cold_indices

