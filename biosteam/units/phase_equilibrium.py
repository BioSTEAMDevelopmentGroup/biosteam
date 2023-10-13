# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains abstract classes for modeling separations in unit operations.

"""
from warnings import warn
from numba import njit
import thermosteam as tmo
from thermosteam import separations as sep
import biosteam as bst
import flexsolve as flx
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from .. import Unit

__all__ = (
    'MultiStageEquilibrium',
)

# %% Equilibrium objects.

@njit(cache=True)
def _vle_phi_K(vapor, liquid):
    F_vapor = vapor.sum()
    F_liquid = liquid.sum()
    phi = F_vapor / (F_vapor + F_liquid)
    y = vapor / F_vapor
    x = liquid / F_liquid
    return phi, y / x 

def _get_specification(name, value):
    if name == 'Duty':
        B = None
        Q = value
    elif name == 'Reflux':
        B = 1 / value
        Q = None
    elif name == 'Boilup':
        B = value
        Q = None
    else:
        raise RuntimeError(f"specification '{name}' not implemented for stage")
    return B, Q


class StageEquilibrium(Unit):
    _N_ins = 0
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    auxiliary_unit_names = ('partition', 'mixer', 'splitters')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                     phases=None, partition_data=None,
                     top_split, bottom_split):
        self._N_outs = 2 + int(top_split) + int(bottom_split)
        Unit.__init__(self, ID, ins, outs, thermo)
        mixer = self.auxiliary(
            'mixer', bst.Mixer, ins=self.ins, conserve_phases=True,
        )
        partition = self.auxiliary(
            'partition', PhasePartition, ins=mixer-0, phases=phases,
            partition_data=partition_data, 
            outs=(
                bst.Stream(None) if top_split else self.outs[0],
                bst.Stream(None) if bottom_split else self.outs[1],
            ),
        )
        self.top_split = top_split
        self.bottom_split = bottom_split
        self.splitters = []
        if top_split:
            self.auxiliary(
                'splitters', bst.Splitter, 
                partition-0, [self.outs[2], self.outs[0]],
                split=top_split,
            )
        if bottom_split:
            self.auxiliary(
                'splitters', bst.Splitter, 
                partition-1, [self.outs[-1], self.outs[1]],
                split=bottom_split, 
            )
    
    def add_feed(self, stream):
        self.ins.append(stream)
        self.mixer.ins.append(
            self.auxlet(
                stream
            )
        )
        
    def set_specification(self, B, Q):
        partition = self.partition
        partition.B = B
        partition.Q = Q
    
    @property
    def extract(self):
        return self.outs[0]
    @property
    def raffinate(self):
        return self.outs[1]
    @property
    def extract_side_draw(self):
        if self.top_split: return self.outs[2]
    @property
    def raffinate_side_draw(self):
        if self.bottom_split: return self.outs[-1]
    
    @property
    def vapor(self):
        return self.outs[0]
    @property
    def liquid(self):
        return self.outs[1]
    @property
    def vapor_side_draw(self):
        if self.top_split: return self.outs[2]
    @property
    def liquid_side_draw(self):
        if self.bottom_split: return self.outs[-1]
    
    def _run(self):
        self.mixer._run()
        self.partition._run()
        for i in self.splitters: i._run()


class PhasePartition(Unit):
    _N_ins = 1
    _N_outs = 2
    strict_infeasibility_check = False
    
    def _init(self, phases, partition_data):
        self.partition_data = partition_data
        self.phases = phases
        self.solvent = None
        self.phi = None
        self.IDs = None
        self.K = None
        self.B = None
        self.Q = 0.
    
    @property
    def extract(self):
        return self.outs[0]
    @property
    def raffinate(self):
        return self.outs[1]
    
    @property
    def vapor(self):
        return self.outs[0]
    @property
    def liquid(self):
        return self.outs[1]
    
    def _run(self, stacklevel=1, P=None, solvent=None, update=True,
             couple_energy_balance=True):
        if solvent is None: solvent = self.solvent
        else: self.solvent = solvent
        for i, j in zip(self.outs, self.phases): i.phase = j
        if update:
            ms = tmo.MultiStream.from_streams(self.outs)
            ms.copy_like(self.feed)
        else:
            ms = self.feed.copy()
            ms.phases = self.phases
        top, bottom = ms
        partition_data = self.partition_data
        if partition_data:
            self.K = K = partition_data['K']
            self.IDs = IDs = partition_data['IDs']
            args = (IDs, K, self.phi or partition_data['phi'], 
                    partition_data.get('extract_chemicals'),
                    partition_data.get('raffinate_chemicals'),
                    self.strict_infeasibility_check, stacklevel+1)
            self.phi = sep.partition(ms, top, bottom, *args)
            self.T = ms.T
        else:
            if 'g' in ms.phases:
                if couple_energy_balance:
                    B = self.B
                    Q = self.Q
                    if B is None: 
                        H = ms.H + Q
                        V = None
                    else:
                        H = None
                        # B = V / (1 - V)
                        # B(1 - V) = V
                        # B - BV - V = 0
                        # -V(1 + B) + B = 0
                        V = B / (1 + B)
                        phi = V
                    ms.vle(P=P or ms.P, H=H, V=V)
                    index = ms.vle._index
                    IDs = ms.chemicals.IDs
                    IDs = tuple([IDs[i] for i in index])
                    L_mol = ms.imol['l', IDs]
                    L_total = L_mol.sum()
                    if L_total: 
                        x_mol = L_mol / L_total
                    else:
                        x_mol = 1
                    V_mol = ms.imol['g', IDs]
                    V_total = V_mol.sum()
                    if V_total: 
                        y_mol = V_mol / V_total
                    else:
                        y_mol = 0
                    K_new = y_mol / x_mol
                    self.phi = V or ms['g'].F_mol / ms.F_mol
                    self.T = ms.T
                else:
                    top, bottom = self.outs
                    if bottom.isempty():
                        p = top.dew_point_at_P(P, IDs=self.IDs)
                    else:
                        p = bottom.bubble_point_at_P(P, IDs=self.IDs)
                    # TODO: Note that solution decomposition method is bubble point
                    x = p.x
                    x[x == 0] = 1.
                    K_new = p.y / p.x
                    self.T = p.T
                    IDs = p.IDs
                    # DO NOT DELETE: Possibly another decomposition method
                    # B = self.B
                    # Q = self.Q
                    # if B is None: 
                    #     H = ms.H + Q
                    # else:
                    #     H = None
                    # ms.vle(P=P or ms.P, H=H, V=B)
                    # IDs = [i.ID for i in ms.vle_chemicals]
                    # fg = ms.imol['g', IDs]
                    # Fg = fg.sum()
                    # fl = ms.imol['l', IDs]
                    # Fl = fl.sum()
                    # self.K = (fg / Fg if Fg else 0.) / (fl / Fl if Fl else 1e-16)
                    # self.phi = Fg / (Fg + Fl)
                    # for i in self.outs: i.T = ms.T
            else:
                eq = ms.lle
                lle_chemicals, K_new, phi = eq(T=ms.T, P=P, top_chemical=solvent, update=update)
                self.phi = phi
                self.T = ms.T
                IDs = tuple([i.ID for i in lle_chemicals])
            IDs_last = self.IDs
            if IDs_last and IDs_last != IDs:
                Ks = self.K
                for ID, K in zip(IDs, K_new): Ks[IDs_last.index(ID)] = K
            else:
                self.K = K_new
                self.IDs = IDs
            

    
class MultiStageEquilibrium(Unit):
    """
    Create a MultiStageEquilibrium object that models counter-current 
    equilibrium stages.
    
    Parameters
    ----------
    N_stages : int
        Number of stages.
    feed_stages : tuple[int]
        Respective stage where feeds enter. Defaults to (0, -1).
    partition_data : {'IDs': tuple[str], 'K': 1d array}, optional
        IDs of chemicals in equilibrium and partition coefficients (molar 
        composition ratio of the extract over the raffinate or vapor over liquid). If given,
        The mixer-settlers will be modeled with these constants. Otherwise,
        partition coefficients are computed based on temperature and composition.
    solvent : str
        Name of main chemical in the solvent.
        
    Examples
    --------
    Simulate 2-stage extraction of methanol from water using octanol:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=500, Methanol=50)
    >>> solvent = bst.Stream('solvent', Octanol=500)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, solvent], phases=('L', 'l'))
    >>> MSE.simulate()
    >>> MSE.extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.83
    >>> MSE.extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    0.99
    >>> MSE.raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    0.82
    
    Simulate 10-stage extraction with user defined partition coefficients:
    
    >>> import numpy as np
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'])
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=10, ins=[feed, solvent], phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.451e-01, 1.380e+00, 2.958e+03]),
    ...         'IDs': ('Water', 'Methanol', 'Octanol'),
    ...         'phi': 0.5899728891780545, # Initial phase fraction guess. This is optional.
    ...     }
    ... )
    >>> MSE.simulate()
    >>> MSE.extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.99
    >>> MSE.extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    0.99
    >>> MSE.raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    0.82
    
    Because octanol and water do not mix well, it may be a good idea to assume
    that these solvents do not mix at all:
        
    >>> MSE = bst.MultiStageEquilibrium(N_stages=20, ins=[feed, solvent], phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSE.simulate()
    >>> MSE.extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.99
    >>> MSE.extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    1.0
    >>> MSE.raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    1.0
       
    Simulate with a feed at the 4th stage:
    
    >>> dilute_feed = bst.Stream('dilute_feed', Water=100, Methanol=2)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed, dilute_feed, solvent], 
    ...     feed_stages=[0, 3, -1],
    ...     phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSE.simulate()
    >>> MSE.extract.imol['Methanol'] / (feed.imol['Methanol'] + dilute_feed.imol['Methanol']) # Recovery
    1.0
    
    Simulate with a 60% extract side draw at the 2nd stage:
    
    >>> MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed, solvent],                         
    ...     top_side_draws={1: 0.6},
    ...     phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSE.simulate()
    >>> extract, raffinate, extract_side_draw, *raffinate_side_draws = MSE.outs
    >>> (extract.imol['Methanol'] + extract_side_draw.imol['Methanol']) / feed.imol['Methanol'] # Recovery
    1.0
    
    >>> Simulate stripping column with 2 stages
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['AceticAcid', 'EthylAcetate', 'Water', 'MTBE'], cache=True)
    >>> feed = bst.Stream('feed', Water=75, AceticAcid=5, MTBE=20, T=320)
    >>> steam = bst.Stream('steam', Water=100, phase='g', T=390)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, steam], feed_stages=[0, -1],
    ...     outs=['vapor', 'liquid'],
    ...     phases=('g', 'l'),
    ... )
    >>> MSE.simulate()
    >>> MSE.vapor.imol['MTBE'] / feed.imol['MTBE']
    1.0
    >>> MSE.vapor.imol['Water'] / feed.imol['Water']
    0.579
    >>> MSE.vapor.imol['AceticAcid'] / feed.imol['AceticAcid']
    0.343
    
    # This feature is not yet ready for users
    # Simulate distillation column with 9 stages, a 0.673 reflux ratio, 
    # 2.57 boilup ratio, and feed at stage 4:
    # >>> import biosteam as bst
    # >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    # >>> feed = bst.Stream('feed', Ethanol=80, Water=100, T=80.215 + 273.15)
    # >>> MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed], feed_stages=[3],
    # ...     outs=['vapor', 'liquid'],
    # ...     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.58)},
    # ...     phases=('g', 'l'),
    # ...     couple_energy_balance=False,
    # ... )
    # >>> MSE.simulate()
    
    """
    _N_ins = 2
    _N_outs = 2
    default_maxiter = 20
    default_molar_tolerance = 0.1
    default_relative_molar_tolerance = 0.001
    auxiliary_unit_names = (
        'stages',
    )
    def __init__(self,  ID='', ins=None, outs=(), thermo=None, 
                 feed_stages=None, top_side_draws=None, bottom_side_draws=None,
                 **kwargs):
        if top_side_draws is None: top_side_draws = {}
        elif not isinstance(top_side_draws, dict): top_side_draws = dict(top_side_draws)
        if bottom_side_draws is None: bottom_side_draws = {}
        elif not isinstance(bottom_side_draws, dict): bottom_side_draws = dict(bottom_side_draws)
        if feed_stages is None: feed_stages = (0, -1)
        self._N_ins = len(feed_stages)
        self._N_outs = 2 + len(top_side_draws) + len(bottom_side_draws)
        Unit.__init__(self, ID, ins, outs, thermo, 
                      feed_stages=feed_stages,
                      top_side_draws=top_side_draws, 
                      bottom_side_draws=bottom_side_draws,
                      **kwargs)
    
    def _init(self,
              N_stages, feed_stages, top_side_draws, bottom_side_draws,
              phases=None, P=101325, stage_specifications=None, 
              partition_data=None, solvent=None, use_cache=None):
        # For VLE look for best published algorithm (don't try simple methods that fail often)
        if phases is None: phases = ('g', 'l')
        if stage_specifications is None: stage_specifications = {}
        elif not isinstance(stage_specifications, dict): stage_specifications = dict(stage_specifications)
        self.multi_stream = tmo.MultiStream(None, P=P, phases=phases, thermo=self.thermo)
        self.N_stages = N_stages
        self.P = P
        phases = self.multi_stream.phases # Corrected order
        top_mark = 2 + len(top_side_draws)
        tsd_iter = iter(self.outs[2:top_mark])
        bsd_iter = iter(self.outs[top_mark:])
        last_stage = None
        self._asplit = asplits = -np.ones(N_stages)
        self._bsplit = bsplits = asplits.copy()
        self.stages = stages = []
        for i in range(N_stages):
            if last_stage is None:
                feed = ()
            else:
                feed = last_stage-1
            outs = []
            if i == 0:
                outs.append(
                    self-0, # extract or vapor
                )
            else:
                outs.append(bst.Stream(None))
            if i == N_stages - 1: 
                outs.append(
                    self-1 # raffinate or liquid
                )
            else:
                outs.append(
                    None
                )
            if i in top_side_draws:
                outs.append(next(tsd_iter))
                top_split = top_side_draws[i]
                asplits[i] += top_split 
            else: 
                top_split = 0
            if i in bottom_side_draws:
                outs.append(next(bsd_iter))
                bottom_split = bottom_side_draws[i]
                bsplits[i] += bottom_split
            else: 
                bottom_split = 0
            
            new_stage = self.auxiliary(
                'stages', StageEquilibrium, phases=phases,
                ins=feed,
                outs=outs,
                partition_data=partition_data,
                top_split=top_split,
                bottom_split=bottom_split,
            )
            if last_stage:
                last_stage.add_feed(new_stage-0)
            last_stage = new_stage
        for feed, stage in zip(self.ins, feed_stages):
            stages[stage].add_feed(self.auxlet(feed))
        self.solvent_ID = solvent
        self.partition_data = partition_data
        self.feed_stages = feed_stages
        self.top_side_draws = top_side_draws
        self.bottom_side_draws = bottom_side_draws
        
        #: dict[int, tuple(str, float)] Specifications for VLE by stage
        self.stage_specifications = stage_specifications
        for i, (name, value) in stage_specifications.items():
            B, Q = _get_specification(name, value)
            stages[i].set_specification(B=B, Q=Q)
            
        #: [int] Maximum number of iterations.
        self.maxiter = self.default_maxiter

        #: [float] Molar tolerance (kmol/hr)
        self.molar_tolerance = self.default_molar_tolerance

        #: [float] Relative molar tolerance
        self.relative_molar_tolerance = self.default_relative_molar_tolerance
        
        self.use_cache = True if use_cache else False
                
    
    @property
    def extract(self):
        return self.outs[0]
    @property
    def raffinate(self):
        return self.outs[1]
    @property
    def extract_side_draws(self):
        return self.top_side_draws
    @property
    def raffinate_side_draws(self):
        return self.bottom_side_draws
    @property
    def vapor(self):
        return self.outs[0]
    @property
    def liquid(self):
        return self.outs[1]
    @property
    def vapor_side_draws(self):
        return self.top_side_draws
    @property
    def liquid_side_draws(self):
        return self.bottom_side_draws
    
    def correct_overall_mass_balance(self):
        outmol = sum([i.mol for i in self.outs])
        inmol = sum([i.mol for i in self.ins])
        factor = inmol / outmol
        for i in self.outs: i.mol *= factor
    
    def material_errors(self):
        errors = []
        stages = self.stages
        IDs = self.multi_stream.chemicals.IDs
        for stage in stages:
            errors.append(
                sum([i.imol[IDs] for i in stage.ins]) - sum([i.imol[IDs] for i in stage.outs])
            )
        return pd.DataFrame(errors, columns=IDs)
    
    def update(self, top_flow_rates):
        top, bottom = self.multi_stream.phases
        flow_tol = -1e-6 * self.multi_stream.mol
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        index = self._update_index
        top_flow_rates[top_flow_rates < 0.] = 0.
        has_infeasible_flow = True
        while has_infeasible_flow:
            has_infeasible_flow = False
            for i in range_stages:
                stage = stages[i]
                partition = stage.partition
                s_top, _ = partition.outs
                s_top.mol[index] = top_flow_rates[i]
                if stage.top_split: stage.splitters[0]._run()
            for i in range_stages:
                stage = stages[i]
                partition = stage.partition
                s_top, s_bottom = partition.outs
                bottom_flow = sum([i.mol for i in stage.ins]) - s_top.mol
                mask = bottom_flow < 0.
                if mask.any():
                    has_infeasible_flow = (bottom_flow[mask] < flow_tol[mask]).any()
                    if has_infeasible_flow:
                        infeasible_index, = np.where(mask[index])
                        # TODO: Find algebraic solution to keeping top flow rates within feasible region.
                        # This is only a temporary solution.
                        infeasible_flow = bottom_flow[mask]
                        top_flow_rates[i, infeasible_index] += infeasible_flow
                        break
                    else:
                        bottom_flow[mask] = 0.
                s_bottom.mol[:] = bottom_flow
                if stage.bottom_split: stage.splitters[-1]._run()
        self.correct_overall_mass_balance()
            
    def _run(self):
        f = self.multi_stage_equilibrium_iter
        top_flow_rates = self.hot_start()
        top_flow_rates = flx.conditional_wegstein(f, top_flow_rates)
        self.update(top_flow_rates)
    
    def hot_start(self):
        self.iter = 1
        ms = self.multi_stream
        feeds = self.ins
        feed_stages = self.feed_stages
        stages = self.stages
        partitions = [i.partition for i in stages]
        N_stages = self.N_stages
        top_phase, bottom_phase = ms.phases
        eq = 'vle' if top_phase == 'g' else 'lle'
        ms.mix_from(feeds)
        ms.P = self.P
        if eq == 'lle':
            self.solvent_ID = solvent_ID = self.solvent_ID or feeds[-1].main_chemical
        data = self.partition_data
        if data:
            top_chemicals = data.get('extract_chemicals') or data.get('vapor_chemicals')
            bottom_chemicals = data.get('raffinate_chemicals') or data.get('liquid_chemicals')
        if self.use_cache and all([not i.multi_stream.isempty() for i in partitions]): # Use last set of data
            if eq == 'lle':
                IDs = data['IDs'] if data else [i.ID for i in ms.lle_chemicals]
                for i in partitions: 
                    if i.IDs != IDs:
                        i.IDs = IDs
                        i._run(P=self.P, solvent=solvent_ID, update=False,
                               couple_energy_balance=False)
                phase_ratios = np.array(
                    [i.phi / (1 - i.phi) for i in partitions],
                    dtype=float
                )
            else:
                IDs = data['IDs'] if data else [i.ID for i in ms.vle_chemicals]
                for i in partitions: 
                    if i.IDs != IDs:
                        i.IDs = IDs
                        i._run(P=self.P, update=False, 
                               couple_energy_balance=False)
                phase_ratios = self.get_vle_phase_ratios()
            IDs = tuple(IDs)
            index = ms.chemicals.get_index(IDs)
            partition_coefficients = np.array([i.K for i in partitions], dtype=float)
            N_chemicals = partition_coefficients.shape[1]
        else:
            if data: 
                top, bottom = ms
                IDs = data['IDs']
                K = data['K']
                phi = data.get('phi') or top.imol[IDs].sum() / ms.imol[IDs].sum()
                data['phi'] = phi = sep.partition(ms, top, bottom, IDs, K, phi,
                                                  top_chemicals, bottom_chemicals)
                index = ms.chemicals.get_index(IDs)
            elif eq == 'lle':
                lle = ms.lle
                lle(ms.T, top_chemical=solvent_ID)
                IDs = tuple([i.ID for i in lle._lle_chemicals])
                index = ms.chemicals.get_index(IDs)
                K = lle._K
                phi = lle._phi
            else:
                P = self.P
                # TODO: Set up much better initial conditions
                if self.stage_specifications:
                    IDs = tuple([i.ID for i in ms.vle_chemicals])
                    dp = ms.dew_point_at_P(P=P, IDs=IDs)
                    T_top = dp.T
                    bp = ms.bubble_point_at_P(P=P, IDs=IDs)
                    T_bot = bp.T
                    dT_stage = (T_bot - T_top) / N_stages
                    for i in range(N_stages):
                        for j in stages[i].partition.outs:
                            j.T = T_top - i * dT_stage
                    index = ms.chemicals.get_index(IDs)
                    phi = 0.5
                    K = bp.y / bp.z
                else:
                    vle = ms.vle
                    vle(H=ms.H, P=P)
                    index = vle._index
                    IDs = ms.chemicals.IDs
                    IDs = tuple([IDs[i] for i in index])
                    L_mol = ms.imol['l', IDs]
                    x_mol = L_mol / L_mol.sum()
                    V_mol = ms.imol['g', IDs]
                    y_mol = V_mol / V_mol.sum()
                    K = y_mol / x_mol
                    phi = ms.V
                    T = ms.T
                    for i in range(N_stages):
                        for i in stages[i].partition.outs:
                            i.T = T
            for i in partitions: 
                i.IDs = IDs
                i.K = K
                i.phi = phi
            N_chemicals = len(index)
            phase_ratios = np.ones(N_stages) * (phi / (1 - phi))
            partition_coefficients = np.ones([N_stages, N_chemicals]) * K[np.newaxis, :]
        if data:
            if top_chemicals:
                top_side_draws = self.top_side_draws
                F = np.zeros([N_stages, len(top_chemicals)])
                top_flow_rates = F.copy()
                for feed, stage in zip(feeds, feed_stages):
                    F[stage] = feed.imol[top_chemicals]
                A = np.eye(N_stages)
                for j, ID in enumerate(top_chemicals):
                    Aj = A.copy()
                    f = F[:, j]
                    for i in range(N_stages - 1):
                        Aj[i, i+1] = -1 
                    for i, value in top_side_draws.items():
                        Aj[i-1, i] *= (1 - value)    
                    top_flow_rates[:, j] = np.linalg.solve(Aj, f)
                for partition, a in zip(partitions, top_flow_rates):
                    partition.outs[0].imol[top_chemicals] = a
                for i in top_side_draws:
                    for s in stages[i].splitters: s._run()
            if bottom_chemicals:
                bottom_side_draws = self.bottom_side_draws
                F = np.zeros([N_stages, len(bottom_chemicals)])
                bottom_flow_rates = F.copy()
                for feed, stage in zip(feeds, feed_stages):
                    F[stage] = feed.imol[bottom_chemicals]
                A = np.eye(N_stages)
                for j, ID in enumerate(bottom_chemicals):
                    Aj = A.copy()
                    f = F[:, j]
                    for i in range(1, N_stages):
                        Aj[i, i-1] = -1 
                    for i, value in bottom_side_draws.items():
                        Aj[i+1, i] *= (1 - value)    
                    bottom_flow_rates[:, j] = np.linalg.solve(Aj, f)
                for partition, b in zip(partitions, bottom_flow_rates):
                    partition.outs[1].imol[bottom_chemicals] = b
                for i in bottom_side_draws:
                    for s in stages[i].splitters: s._run()
        self.feed_flows = feed_flows = np.zeros([N_stages, N_chemicals])
        for feed, stage in zip(feeds, feed_stages):
            feed_flows[stage, :] += feed.mol[index]
        top_flow_rates = flow_rates_for_multi_stage_equilibrium(
            phase_ratios, partition_coefficients, feed_flows, self._asplit, self._bsplit,
        )
        self._iter_args = (feed_flows, self._asplit, self._bsplit)
        self._update_index = index
        self._N_chemicals = N_chemicals
        return top_flow_rates
    
    def get_vle_phase_ratios(self):
        # ENERGY BALANCE
        # Intermediate stage:
        # Hv1*V1 + Hl1*L1 - Hv2*V2 - Hl0*L0 = Q1
        # Hv1*L1*B1 + Hl1*L1 - Hv2*L2*B2 - Hl0*L0 = Q1
        # Hv1*L1*B1 - Hv2*L2*B2 = Q1 + Hl0*L0 - Hl1*L1
        # B1 = (Q1 + Hl0*L0 - Hl1*L1 + Hv2*L2*B2) / (Hv1*L1)
        # Top stage:
        # Hv1*V1 + Hl1*L1 - Hv2*V2 = Q1
        # Hv1*L1*B1 + Hl1*L1 - Hv2*L2*B2 = Q1
        # Hv1*L1*B1 + Hl1*L1 - Hv2*L2*B2 = Q1 
        # B1 = (Q1 - Hl1*L1 + Hv2*L2*B2) / (Hv1*L1)
        # Bottom stage:
        # Hv1*V1 + Hl1*L1 - Hl0*L0 = Q1
        # Hv1*L1*B1 + Hl1*L1 - Hl0*L0 = Q1
        # Hv1*L1*B1 = Q1 + Hl0*L0 - Hl1*L1
        # B1 = Q1 + Hl0*L0 - Hl1*L1 / (Hv1*L1)
        # If top stage (or bottom stage) B1 given, do not include include in equations.
        # Second to last (or fist) stage equations do not change. 
        stages = self.stages
        N_stages = self.N_stages
        phase_ratios = np.zeros(N_stages)
        last_stage = vapor_last = liquid_last = B_last = None
        for i in range(N_stages - 1, -1, -1):
            # breakpoint()
            stage = stages[i]
            partition = stage.partition
            B = partition.B
            if B is not None:
                phase_ratios[i] = B_last = B
                continue
            Q = partition.Q or 0
            vapor, liquid, *_ = stage.outs
            Q -= liquid.H
            if last_stage:
                Q += vapor_last.h * liquid_last.F_mol * B_last
            if i != 0:
                next_stage = stages[i - 1] # Going up
                vapor_next, liquid_next, *_ = next_stage.outs
                Q += liquid_next.H
            adjacent_stages = [i for i in (next_stage, last_stage) if i]
            other_feeds = [i for i in stage.ins if i.source not in adjacent_stages]
            Q += sum([i.H for i in other_feeds])
            B = Q / (vapor.h * liquid.F_mol)
            vapor_last = vapor
            liquid_last = liquid
            phase_ratios[i] = B_last = B
            last_stage = stage
        return phase_ratios
    
    def multi_stage_equilibrium_iter(self, top_flow_rates):
        self.iter += 1
        self.update(top_flow_rates)
        stages = self.stages
        P = self.P
        eq = 'vle' if self.multi_stream.phases[0] == 'g' else 'lle'
        if eq == 'vle': 
            for n, i in enumerate(stages): 
                i.mixer._run()
                i.partition._run(P=self.P, update=False, 
                                 couple_energy_balance=False)
            for i, j in enumerate(self.get_vle_phase_ratios()):
                stages[i].partition.phi = 1 / (1 + j)
        else:
            for i in stages: 
                i.mixer._run()
                i.partition._run(P=P, solvent=self.solvent_ID, update=False, 
                                 couple_energy_balance=False)
        phase_ratios = []
        partition_coefficients = []
        Ts = []
        almost_one = 1. - 1e-16
        almost_zero = 1e-16
        N_stages = self.N_stages
        index = []
        for i in range(N_stages):
            partition = stages[i].partition
            phi = partition.phi
            if almost_zero < phi < almost_one:
                index.append(i)
                phase_ratios.append(phi / (1 - phi))
                partition_coefficients.append(partition.K)
                Ts.append(partition.T)
        N_ok = len(index)
        if len(index) == N_stages:
            phase_ratios = np.array(phase_ratios)
            partition_coefficients = np.array(partition_coefficients)
            Ts = np.array(Ts)
        elif N_ok > 1:
            all_index = np.arange(N_stages)
            spline = UnivariateSpline(index, phase_ratios, k=1, ext='const')
            phase_ratios = spline(all_index)
            spline = UnivariateSpline(index, Ts, k=1, ext='const')
            Ts = spline(all_index)
            N_chemicals = self._N_chemicals
            all_partition_coefficients = np.zeros([N_stages, N_chemicals])
            for i in range(N_chemicals):
                spline = UnivariateSpline(index, [stage[i] for stage in partition_coefficients], k=1, ext='const')
                all_partition_coefficients[:, i] = spline(all_index)
            partition_coefficients = all_partition_coefficients
        elif N_ok == 1:
            phase_ratios = np.array(N_stages * phase_ratios)
            partition_coefficients = np.array(N_stages * partition_coefficients)
            Ts = np.array(N_stages * Ts)
        for T, stage in zip(Ts, stages): 
            for i in stage.partition.outs: i.T = T
        new_top_flow_rates = flow_rates_for_multi_stage_equilibrium(
            phase_ratios, partition_coefficients, *self._iter_args,
        )
        mol = top_flow_rates[0] 
        mol_new = new_top_flow_rates[0]
        mol_errors = abs(mol - mol_new)
        if mol_errors.any():
            mol_error = mol_errors.max()
            if mol_error > 1e-12:
                nonzero_index, = (mol_errors > 1e-12).nonzero()
                mol_errors = mol_errors[nonzero_index]
                max_errors = np.maximum.reduce([abs(mol[nonzero_index]), abs(mol_new[nonzero_index])])
                rmol_error = (mol_errors / max_errors).max()
                not_converged = (
                    self.iter < self.maxiter and (mol_error > self.molar_tolerance
                     or rmol_error > self.relative_molar_tolerance)
                )
            else:
                not_converged = False
        else:
            not_converged = False
        return new_top_flow_rates, not_converged


# %% General functional algorithms based on MESH equations to solve multi-stage 

@njit(cache=True)
def solve_TDMA(a, b, c, d): # Tridiagonal matrix solver
    """
    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    """
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        m = a[i] / b[i]
        b[inext] = b[inext] - m * c[i] 
        d[inext] = d[inext] - m * d[i]
        
    b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        b[i] = (d[i] - c[i] * b[i+1]) / b[i]
    
    return b

@njit(cache=True)
def flow_rates_for_multi_stage_equilibrium(
        phase_ratios,
        partition_coefficients, 
        feed_flows,
        asplit,
        bsplit,
    ):
    """
    Solve b-phase flow rates for a single component across equilibrium stages with side draws. 

    Parameters
    ----------
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a (extract or vapor) and b (raffinate or liquid) leaving the stage 
        respectively.
    partition_coefficients : Iterable[1d array]
        Partition coefficients with stages by row and components by column.
        The partition coefficient for a component in a given stage is defined 
        as x_a / x_b; where x_a and x_b are the fraction of the component in 
        phase a (extract or vapor) and b (raffinate or liquid) leaving the stage.
    feed_flows : Iterable [1d array]
        Flow rates of all components feed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    phase_ratios[phase_ratios < 1e-16] = 1e-16
    phase_ratios[phase_ratios > 1e16] = 1e16
    N_stages, N_chemicals = partition_coefficients.shape
    phase_ratios = np.expand_dims(phase_ratios, -1)
    component_ratios = 1 / (phase_ratios * partition_coefficients)
    b = 1. +  component_ratios
    a = b.copy()
    c = b.copy()
    d = feed_flows.copy()
    for i in range(N_stages-1):
        c[i] = bsplit[i + 1]
        a[i] = asplit[i] *  component_ratios[i] 
    return solve_TDMA(a, b, c, d)