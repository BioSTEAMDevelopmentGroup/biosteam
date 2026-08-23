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
import heapq
import numpy as np
import biosteam as bst
from warnings import warn

__all__ = ('StreamLifeCycle', 'ProblemTable', 'problem_table',
           'synthesize_network', 'plot_pinch_diagram')

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

def _stream_H_at_boundaries(stream_in, H_in, H_out, T_lo, T_hi, Ts, shift,
                            stream_label):
    """
    Enthalpies [kJ/hr] of one monotone stream at the grid boundaries
    `Ts` (shifted scale, descending, all within [T_lo, T_hi]).

    Exact at the stream's own end points (H_in/H_out as given) by
    *position*: `Ts[0]` and `Ts[-1]` are the stream's own T_hi/T_lo (every
    monotone stream has `T_hi > T_lo` strictly, so `Ts` always has at least
    these two entries) and are assigned H_in/H_out directly, without a float
    comparison. In between, the inlet copy is flashed at the *real*
    temperature `T + shift` and the result is clipped to
    [min(H_in, H_out), max(H_in, H_out)] so that a non-equilibrium outlet
    (e.g. a column reboiler/condenser product) can never inflate an
    interval. A single copy is walked down the grid so each VLE is
    warm-started from the previous boundary; `stream_label` (the inlet
    stream's own ID) identifies the stream in the VLE-failure warning.
    """
    assert Ts.size >= 2, (
        "boundary grid for a monotone stream must include both its own "
        "end points"
    )
    H_lo, H_hi = sorted((H_in, H_out))
    H_top, H_bottom = (H_in, H_out) if H_in > H_out else (H_out, H_in)
    Hs = np.empty(Ts.size)
    Hs[0] = H_top
    Hs[-1] = H_bottom
    stream = stream_in.copy()
    for k in range(1, Ts.size - 1):
        T = Ts[k]
        T_real = T + shift
        try:
            stream.vle(T=T_real, P=stream.P)
            H = stream.H
        except Exception as error:
            warn(f"could not solve VLE for stream {stream_label!r} at "
                 f"{T_real:.2f} K ({error!r}); interpolating enthalpy "
                 "linearly in temperature for the problem table",
                 RuntimeWarning)
            # restart the warm start from a clean copy so the failed flash
            # does not leave `stream` in a bad state for the next boundary
            stream = stream_in.copy()
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
        hot, - for cold), the cascade `residual` (n) *leaving* each
        boundary (i.e. after its point loads), `hot_util_load`,
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
    is residual[k] = sum(point_H[:, :k+1]) + sum(interval_H[:, :k]), the
    heat *leaving* boundary Ts[k]. Feasibility must also hold for the heat
    *arriving* at Ts[k] before its point loads are applied,
    arriving[k] = residual[k] - sum(point_H[:, k]), because a source at
    Ts[k] cannot serve a sink above Ts[k]. The minimum over both flows,
    min(residual, arriving), fixes the hot utility target,
    `residual[-1] + hot_util_load` the cold one, and its location the
    pinch. With the per-stream identity above, hot_util_load -
    cold_util_load equals the net heating demand.

    Examples
    --------
    A threshold problem: 1000 kmol/hr of water cooled 400 -> 300 K supplies
    every interval of 900 kmol/hr of water heated 300 -> 390 K, so no hot
    utility is needed and the surplus leaves as cold utility.

    >>> import biosteam as bst
    >>> from biosteam.facilities.hxn.hxn_synthesis import problem_table
    >>> bst.settings.set_thermo(['Water'])
    >>> hot_in = bst.Stream(Water=1000., T=400., P=5e5, phase='l', units='kmol/hr')
    >>> hot_out = hot_in.copy(); hot_out.vle(T=300., P=5e5)
    >>> cold_in = bst.Stream(Water=900., T=300., P=5e5, phase='l', units='kmol/hr')
    >>> cold_out = cold_in.copy(); cold_out.vle(T=390., P=5e5)
    >>> table = problem_table([hot_in, cold_in], [hot_out, cold_out],
    ...                       [True, False], 5.)
    >>> table.Ts
    array([395., 390., 300., 295.])
    >>> round(table.hot_util_load, 3)
    0.0
    >>> round(table.cold_util_load, -1)
    1445550.0
    >>> table.pinch_T
    395.0
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
                                         T_lo[j], T_hi[j], Ts[idx], shift[j],
                                         streams_inlet[j].ID)
            interval_H[j, idx[:-1]] = sign[j] * (Hs[:-1] - Hs[1:])
        else:
            k = np.searchsorted(-Ts, -T_hi[j])
            point_H[j, k] = sign[j] * abs(H_out[j] - H_in[j])
    point_total = point_H.sum(axis=0)
    residual = np.cumsum(
        point_total + np.concatenate([[0.], interval_H.sum(axis=0)])
    )
    # heat arriving at each boundary, before that boundary's point loads:
    # a point source at Ts[k] cannot serve sinks above Ts[k], so the cascade
    # must be non-negative both before and after the point loads
    arriving = residual - point_total
    flow = np.minimum(residual, arriving)
    k_pinch = int(np.argmin(flow))
    scale = np.abs(H_out - H_in).sum()
    if -flow[k_pinch] <= 1e-9 * scale:  # threshold problem: no hot utility
        hot_util_load = 0.
        k_pinch = 0
    else:
        hot_util_load = -flow[k_pinch]
    cold_util_load = residual[-1] + hot_util_load
    if cold_util_load < 0.:
        # only reachable in the threshold branch, by at most 1e-9 * scale:
        # absorb the rounding into the hot utility so that
        # hot_util_load - cold_util_load == sum(unit_duty) stays exact
        hot_util_load -= cold_util_load
        cold_util_load = 0.
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
    # Per-stream pinch temperature: where each stream is split between the
    # hot-side and cold-side network designs. A stream already entirely on
    # one side of the process pinch (T_in past pinch_cold_stream_T for a
    # cold stream, or past pinch_hot_stream_T for a hot stream) is not
    # split; its pinch_T is its own T_in. This clause also catches
    # non-monotone streams (T_out on the wrong side of T_in for their duty,
    # e.g. a cold stream whose VLE outlet ends up cooler than it entered):
    # rather than split their problem_table point-load duty across the
    # cascade, they get pinch_T = T_in too, so load_duties assigns their
    # whole duty to a single side (Q_hot_side for a cold stream,
    # Q_cold_side for a hot one).
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



# Pinch diagram

def _order_exchanger_columns(hxs, stream_life_cycles):
    """
    Order heat exchangers left to right so that every stream meets its
    exchangers in flow direction (cold streams flow left to right, hot streams
    right to left). The per-stream stage orders define a precedence graph;
    a topological sort (Kahn's algorithm, ties broken by the given order)
    yields a consistent layout. Contradictory constraints, which would need a
    stream to flow backwards, fall back to the given order.
    """
    hxs = list(hxs)
    position = {hx: i for i, hx in enumerate(hxs)}
    successors = {hx: [] for hx in hxs}
    N_predecessors = {hx: 0 for hx in hxs}
    for life_cycle in stream_life_cycles:
        stages = [i.unit for i in life_cycle.life_cycle if i.unit in position]
        if not life_cycle.cold: stages.reverse()
        for a, b in zip(stages, stages[1:]):
            if b not in successors[a]:
                successors[a].append(b)
                N_predecessors[b] += 1
    ready = [position[hx] for hx in hxs if not N_predecessors[hx]]
    heapq.heapify(ready)
    ordered = []
    while ready:
        hx = hxs[heapq.heappop(ready)]
        ordered.append(hx)
        for other in successors[hx]:
            N_predecessors[other] -= 1
            if not N_predecessors[other]: heapq.heappush(ready, position[other])
    return ordered if len(ordered) == len(hxs) else hxs

def _format_H(H):
    mantissa, exponent = f'{H:.2e}'.split('e')
    return f'{mantissa}E{int(exponent)}'

def _auxiliary_name(unit):
    """
    Return the (dotted) name of an auxiliary unit within its owner, e.g.
    'condenser' or 'evaporators[0].heat_exchanger', or None if the unit is
    not auxiliary.
    """
    owner = unit.owner
    if owner is unit: return None
    def search(parent, prefix):
        for name, aux in parent.get_auxiliary_units_with_names():
            if aux is unit: return prefix + name
            if hasattr(aux, 'get_auxiliary_units_with_names'):
                found = search(aux, prefix + name + '.')
                if found: return found
    return search(owner, '') or unit.ID.lstrip('.')

def _stream_label(hx, show_units, show_auxiliary_units, show_stream_IDs):
    """
    Label of a stream from its original heat exchanger `hx`:
    '<owner> - <auxiliary name> (<inlet stream ID>)', with each part
    optional.
    """
    parts = []
    if show_units: parts.append(hx.owner.ID)
    if show_auxiliary_units:
        auxname = _auxiliary_name(hx)
        if auxname: parts.append(auxname)
    label = ' - '.join(parts)
    if show_stream_IDs:
        ID = hx.ins[0].ID
        if ID: label = f'{label} ({ID})' if label else ID
    return label

def plot_pinch_diagram(stream_life_cycles, inlet_Ts, outlet_Ts,
                       hot_side_HXs, cold_side_HXs, Qmin=1e-3,
                       original_hxs=None, show_units=True,
                       show_auxiliary_units=True, show_stream_IDs=True,
                       show_legend=True, ax=None, file=None, dpi=300):
    """
    Draw a pinch diagram of a synthesized heat exchanger network: cold
    streams (blue, flowing left to right) above hot streams (red, flowing
    right to left), one vertical connector per process heat exchanger with
    its duty, a dashed pinch line separating the cold-side from the
    hot-side exchangers, and circles marking the utility exchangers that
    bring each stream to its outlet temperature.

    Parameters
    ----------
    stream_life_cycles : list[StreamLifeCycle]
        One per stream, as built by HeatExchangerNetwork.
    inlet_Ts, outlet_Ts : array-like
        Stream inlet and outlet temperatures [K], indexed like the life cycles.
    hot_side_HXs, cold_side_HXs : list[HXprocess]
        Process exchangers above and below the pinch.
    Qmin : float, optional
        Utility exchangers with a duty at or below this [kJ/hr] are not marked.
    original_hxs : list[Unit], optional
        The original heat exchanger of each stream (indexed like the life
        cycles). Required for the stream labels below.
    show_units : bool, optional
        Label each stream with the unit operation that owns its original heat
        exchanger (the main unit for auxiliary exchangers).
    show_auxiliary_units : bool, optional
        Label each stream with the name of its original heat exchanger within
        the main unit (e.g. 'condenser'), if it is an auxiliary unit.
    show_stream_IDs : bool, optional
        Label each stream with the ID of the original heat exchanger's inlet.
    show_legend : bool, optional
        Add a legend of the symbols below the diagram.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if not given.
    file : str, optional
        If given, the figure is saved to this path.
    dpi : int, optional
        Resolution used when saving.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes

    Notes
    -----
    Temperatures are shown in degC and heat flows in kJ/hr at the inlet and
    outlet of each stream. Exchanger columns on each side of the pinch are
    ordered so that each stream meets them in flow direction whenever the
    network allows it. Stream labels read '<unit> - <auxiliary> (<stream>)'
    next to the stream index at the inlet.

    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Glycerol'])
    >>> feed1 = bst.Stream('feed1', flow=(8000, 100, 25))
    >>> feed2 = bst.Stream('feed2', flow=(10000, 1000, 10))
    >>> D1 = bst.ShortcutColumn('D1', ins=feed1,
    ...                     outs=('distillate', 'bottoms_product'),
    ...                     LHK=('Methanol', 'Water'),
    ...                     y_top=0.99, x_bot=0.01, k=2,
    ...                     is_divided=True)
    >>> D1_H1 = bst.HXutility('D1_H1', ins = D1.outs[1], T = 300)
    >>> D1_H2 = bst.HXutility('D1_H2', ins = D1.outs[0], T = 300)
    >>> F1 = bst.Flash('F1', ins=feed2,
    ...                outs=('vapor', 'liquid'), V = 0.9, P = 101325)
    >>> HXN = bst.HeatExchangerNetwork('HXN', T_min_app = 5.)
    >>> sys = bst.System.from_units('sys', units=[D1, D1_H1, D1_H2, F1, HXN])
    >>> sys.simulate()
    >>> fig, ax = HXN.plot_pinch_diagram()
    >>> connectors = [i for i in ax.findobj() if (i.get_gid() or '').startswith('HX:')]
    >>> len(connectors) == len(HXN.new_HXs)
    True
    >>> import matplotlib.pyplot as plt
    >>> plt.close(fig)

    """
    import matplotlib.pyplot as plt
    # Artists carry stable gids ('HX:<ID>', 'Util:<ID>', 'Label:<index>')
    # so the drawing can be checked structurally in tests.
    show_labels = show_units or show_auxiliary_units or show_stream_IDs
    if show_labels and original_hxs is None:
        raise ValueError('original_hxs is required to label streams with '
                         'units, auxiliary units, or stream IDs')
    cold_color, hot_color = '#2e6db4', '#d62728'
    cold_bg, hot_bg = '#e6f0fa', '#fbe9e7'
    process_hxs = set(hot_side_HXs) | set(cold_side_HXs)
    # Stream index and stage of each side of every process exchanger, by identity
    hx_streams = {hx: {} for hx in process_hxs}
    for index, life_cycle in enumerate(stream_life_cycles):
        for stage in life_cycle.life_cycle:
            if stage.unit in hx_streams:
                hx_streams[stage.unit][life_cycle.cold] = (index, stage)
    cold_side_HXs = _order_exchanger_columns(cold_side_HXs, stream_life_cycles)
    hot_side_HXs = _order_exchanger_columns(hot_side_HXs, stream_life_cycles)
    columns = cold_side_HXs + hot_side_HXs
    N_cs = len(cold_side_HXs)
    N_columns = len(columns)
    # x layout: 0 stream ends | 1 cold utilities | 2..N_cs+1 cold side |
    # pinch | N_cs+2..N+1 hot side | N+2 hot utilities | N+3 stream ends
    x_start, x_cold_util = 0., 1.
    x_columns = {hx: 2. + i for i, hx in enumerate(columns)}
    x_pinch = N_cs + 1.5
    x_hot_util = N_columns + 2.
    x_end = N_columns + 3.
    # y layout: cold streams on top, hot streams below, duty labels in between
    cold_streams = [i for i, lc in enumerate(stream_life_cycles) if lc.cold]
    hot_streams = [i for i, lc in enumerate(stream_life_cycles) if not lc.cold]
    N_hot = len(hot_streams)
    N_cold = len(cold_streams)
    gap = 2.5
    y = {}
    for k, i in enumerate(hot_streams): y[i] = N_hot - k
    for k, i in enumerate(cold_streams): y[i] = N_hot + gap + N_cold - k
    y_label = N_hot + (gap + 1.) / 2.
    y_top = N_hot + gap + N_cold + 1.
    y_bottom = 0.
    if ax is None:
        fig, ax = plt.subplots(
            figsize=(max(6., 0.75 * (N_columns + 4) + 3.), 0.4 * y_top + 1.)
        )
    else:
        fig = ax.figure
    # Background and pinch line
    x_min, x_max = x_start - 1.8, x_end + 1.8
    ax.axvspan(x_min, x_pinch, color=cold_bg, lw=0, zorder=0)
    ax.axvspan(x_pinch, x_max, color=hot_bg, lw=0, zorder=0)
    ax.axvline(x_pinch, color='k', ls='--', lw=1, zorder=1)
    ax.text(x_min + 0.2, y_bottom + 0.1, 'Cold side', color=cold_color,
            weight='bold', ha='left', va='bottom')
    ax.text(x_max - 0.2, y_bottom + 0.1, 'Hot side', color=hot_color,
            weight='bold', ha='right', va='bottom')
    # Column headers
    header_kwargs = dict(ha='center', va='bottom', weight='bold', fontsize=8)
    for x_T, x_H in ((x_start - 1.3, x_start - 0.6), (x_end + 0.6, x_end + 1.3)):
        ax.text(x_T, y_top, 'T\n[°C]', **header_kwargs)
        ax.text(x_H, y_top, 'H\n[kJ·h$^{-1}$]', **header_kwargs)
    ax.text(x_start - 0.3, y_label, 'ΔH\n[kJ·h$^{-1}$]',
            ha='right', va='center', weight='bold', fontsize=8)
    # Streams
    value_kwargs = dict(ha='center', va='center', fontsize=8)
    for index, life_cycle in enumerate(stream_life_cycles):
        cold = life_cycle.cold
        color = cold_color if cold else hot_color
        yi = y[index]
        stages = life_cycle.life_cycle # never empty: each stream has a utility stage
        H_in = stages[0].H_in
        H_out = stages[-1].H_out
        T_in = inlet_Ts[index] - 273.15
        T_out = outlet_Ts[index] - 273.15
        # T is the outer column on the left and the inner column on the right
        T_left, H_left, T_right, H_right = (
            (T_in, H_in, T_out, H_out) if cold else (T_out, H_out, T_in, H_in)
        )
        x_in, x_out, sign = (x_start, x_end, 1) if cold else (x_end, x_start, -1)
        ax.annotate('', xy=(x_out, yi), xytext=(x_in, yi),
                    arrowprops=dict(arrowstyle='-|>', color=color, lw=1.2,
                                    shrinkA=0, shrinkB=0), zorder=2)
        ax.text(x_start - 1.3, yi, f'{T_left:.1f}', color=color, **value_kwargs)
        ax.text(x_start - 0.6, yi, _format_H(H_left), color=color, **value_kwargs)
        ax.text(x_end + 0.6, yi, f'{T_right:.1f}', color=color, **value_kwargs)
        ax.text(x_end + 1.3, yi, _format_H(H_right), color=color, **value_kwargs)
        ax.text(x_in + sign * 0.3, yi + 0.12, str(index), color=color,
                ha='center', va='bottom', weight='bold', fontsize=9)
        if show_labels:
            label = _stream_label(original_hxs[index], show_units,
                                  show_auxiliary_units, show_stream_IDs)
            ax.text(x_in + sign * 0.6, yi + 0.12, label, color=color,
                    ha='left' if cold else 'right', va='bottom', fontsize=7,
                    zorder=6, gid=f'Label:{index}',
                    bbox=dict(boxstyle='square,pad=0.15', fc='w', ec='none'))
        # Utility exchangers: a cold stream ends in a hot utility (red), a
        # hot stream in a cold utility (blue)
        x_util, util_color = (x_hot_util, hot_color) if cold else (x_cold_util, cold_color)
        for stage in stages:
            unit = stage.unit
            if unit in process_hxs: continue
            if abs(stage.H_out - stage.H_in) <= Qmin: continue
            ax.plot([x_util], [yi], 'o', mfc='w', mec=util_color, mew=1.2,
                    ms=6, zorder=4, gid='Util:' + unit.ID)
    # Process exchangers
    for hx in columns:
        streams = hx_streams[hx]
        if len(streams) != 2:
            warn(f'{hx.ID} is not in exactly one hot and one cold stream '
                 'life cycle; it is not drawn', RuntimeWarning)
            continue
        (i_cold, stage_cold), (i_hot, stage_hot) = streams[True], streams[False]
        x = x_columns[hx]
        Q = abs(stage_hot.H_in - stage_hot.H_out)
        ax.plot([x, x], [y[i_hot], y[i_cold]], '-o', color='k', mfc='w',
                mew=1.2, ms=6, lw=1.2, zorder=3, gid='HX:' + hx.ID)
        ax.text(x, y_label, _format_H(Q), rotation=90, ha='center',
                va='center', fontsize=8, zorder=5,
                bbox=dict(boxstyle='square,pad=0.25', fc='w', ec='k', lw=0.8))
    if show_legend:
        from matplotlib.lines import Line2D
        handles = [
            Line2D([], [], color=cold_color, lw=1.2, marker='>', markevery=[-1],
                   ms=5, label='Cold stream'),
            Line2D([], [], color=hot_color, lw=1.2, marker='<', markevery=[0],
                   ms=5, label='Hot stream'),
            Line2D([], [], color='k', lw=1.2, marker='o', mfc='w', mew=1.2,
                   ms=6, label='Process heat exchange'),
            Line2D([], [], ls='', marker='o', mfc='w', mec=hot_color, mew=1.2,
                   ms=6, label='Hot utility'),
            Line2D([], [], ls='', marker='o', mfc='w', mec=cold_color, mew=1.2,
                   ms=6, label='Cold utility'),
            Line2D([], [], color='k', ls='--', lw=1, label='Pinch'),
        ]
        ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.),
                  ncol=3, fontsize=7, frameon=False, handlelength=2.5,
                  columnspacing=1.5)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_bottom, y_top + 1.2)
    ax.set_axis_off()
    if file: fig.savefig(file, dpi=dpi, bbox_inches='tight')
    return fig, ax
