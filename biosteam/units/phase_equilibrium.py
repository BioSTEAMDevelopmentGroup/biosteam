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
from numba import njit, objmode
import thermosteam as tmo
from thermosteam import separations as sep
import biosteam as bst
import flexsolve as flx
import numpy as np
import pandas as pd
from scipy.optimize import minimize, differential_evolution
from math import inf
from typing import Callable
from scipy.optimize import root
from ..exceptions import Converged
from .. import Unit

__all__ = (
    'SinglePhaseStage',
    'StageEquilibrium',
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
        if value is None: 
            B = None
        else:
            B = inf if value == 0 else 1 / value
        Q = None
    elif name == 'Boilup':
        B = value
        Q = None
    else:
        raise RuntimeError(f"specification '{name}' not implemented for stage")
    return B, Q
        
class SinglePhaseStage(Unit):
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    
    def _init(self, T=None, phase=None):
        self.T = T
        self.phase = phase
        self.stages = [self]
        
    @property
    def phases(self):
        return (self.phase,)
        
    def _run(self):
        outlet = self.outs[0]
        outlet.mix_from(self.ins, energy_balance=False)
        outlet.phase = self.phase
        outlet.T = self.T

    # %% Decoupled phenomena equation oriented simulation
    
    def _get_decoupled_variable(self, variable):
        if variable == 'T': return self.T
        
    def _update_decoupled_variable(self, variable, value): 
        pass
        
    def _create_mass_balance_equations(self):
        outlet = self.outs[0]
        inlets = self.ins
        fresh_inlets = [i for i in inlets if i.isfeed()]
        process_inlets = [i for i in inlets if not i.isfeed()]
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {outlet: ones}
        for i in process_inlets: eq_overall[i] = minus_ones
        return [
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        ]
    
    def _create_linear_equations(self, variable):
        if variable == 'mol':
            return self._create_mass_balance_equations()
        elif variable == 'mol-LLE':
            return [] # TODO: Create case where mixers are connected
        else:
            return []

# %%

class StageEquilibrium(Unit):
    _N_ins = 0
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    auxiliary_unit_names = ('partition', 'mixer', 'splitters')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
            phases, partition_data=None, top_split=0, bottom_split=0,
            B=None, Q=None, top_chemical=None,
        ):
        self._N_outs = 2 + int(top_split) + int(bottom_split)
        self.phases = phases
        self.stages = [self]
        Unit.__init__(self, ID, ins, outs, thermo)
        mixer = self.auxiliary(
            'mixer', bst.Mixer, ins=self.ins, 
        )
        mixer.outs[0].phases = phases
        partition = self.auxiliary(
            'partition', PhasePartition, ins=mixer-0, phases=phases,
            partition_data=partition_data, top_chemical=top_chemical,
            outs=(
                None if top_split else self.outs[0],
                None if bottom_split else self.outs[1],
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
        self.set_specification(B, Q)
    
    @property
    def Q(self):
        return self.partition.Q
    @Q.setter
    def Q(self, Q):
        self.partition.Q = Q
    
    @property
    def B(self):
        return self.partition.B
    @B.setter
    def B(self, B):
        self.partition.B = B
    
    @property
    def B_specification(self):
        return self.partition.B_specification
    @B_specification.setter
    def B_specification(self, B_specification):
        self.partition.B_specification = B_specification
    
    @property
    def T(self):
        return self.partition.T
    @T.setter
    def T(self, T):
        self.partition.T = T
        for i in self.partition.outs: i.T = T
    
    @property
    def K(self):
        return self.partition.K
    @K.setter
    def K(self, K):
        self.partition.K = K
    
    def _update_auxiliaries(self):
        for i in self.splitters: i.ins[0].mix_from(i.outs, energy_balance=False)
        self.mixer.outs[0].mix_from(self.ins, energy_balance=False)
    
    def add_feed(self, stream):
        self.ins.append(stream)
        self.mixer.ins.append(
            self.auxlet(
                stream
            )
        )
        
    def set_specification(self, B, Q):
        if B is None and Q is None: Q = 0.
        partition = self.partition
        partition.B_specification = partition.B = B
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
    @property
    def top_side_draw(self):
        if self.top_split: return self.outs[2]
    @property
    def bottom_side_draw(self):
        if self.bottom_split: return self.outs[-1]
    
    def _run(self):
        self.mixer._run()
        self.partition._run()
        for i in self.splitters: i._run()
    
    # %% Decoupled phenomena equation oriented simulation
    
    def _create_temperature_departure_equation(self):
        # C1dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - H_out + H_in
        coeff = {self: sum([i.C for i in self.outs])}
        for i in self.ins:
            source = i.source
            if not source or source.phases != ('L', 'l'): continue
            coeff[source] = -i.C
        return coeff, (self.Q or 0.) + self.H_in - self.H_out
    
    def _create_phase_ratio_departure_equation(self):
        # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
        vapor, liquid = self.partition.outs
        coeff = {}
        if vapor.isempty():
            liquid.phase = 'g'
            coeff[self] = liquid.H
            liquid.phase = 'l'
        else:
            coeff[self] = vapor.h * liquid.F_mol
        for i in self.ins:
            if i.phase != 'g': continue
            source = i.source
            if (not source
                or source.phases != ('g', 'l')
                or source.B_specification is not None):
                continue
            vapor, liquid, *_ = source.outs
            if vapor.isempty():
                liquid.phase = 'g'
                coeff[source] = liquid.H
                liquid.phase = 'l'
            else:
                coeff[source] = -vapor.h * liquid.F_mol
        return coeff, (self.Q or 0.) + self.H_in - self.H_out
    
    def _create_lle_mass_balance_equations(self):
        top_split = self.top_split
        bottom_split = self.bottom_split
        inlets = self.ins
        fresh_inlets = []
        process_inlets = []
        for i in inlets:
            if not i.isfeed() and i.source.phases == ('L' 'l'):
                process_inlets.append(i)
            else:
                fresh_inlets.append(i)
        top, bottom, *_ = self.outs
        top_side_draw = self.top_side_draw
        bottom_side_draw = self.bottom_side_draw
        equations = []
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {}
        S = self.K * self.B
        for i in self.outs: eq_overall[i] = ones
        for i in process_inlets: eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        )
        
        # Top to bottom flows
        eq_outs = {}
        eq_outs[top] = ones
        eq_outs[bottom] = -S
        equations.append(
            (eq_outs, zeros)
        )
        
        # Top split flows
        if top_side_draw:
            eq_top_split = {
                top_side_draw: top_split * ones,
                top: top_split + minus_ones,
            }
            equations.append(
                (eq_top_split, zeros)
            )
        
        # Bottom split flows
        if bottom_side_draw:
            eq_bottom_split = {
                bottom_side_draw: bottom_split * ones,
                bottom: bottom_split + minus_ones,
            }
            equations.append(
                (eq_bottom_split, zeros)
            )
        
        return equations
    
    def _create_mass_balance_equations(self):
        top_split = self.top_split
        bottom_split = self.bottom_split
        inlets = self.ins
        fresh_inlets = [i for i in inlets if i.isfeed()]
        process_inlets = [i for i in inlets if not i.isfeed()]
        top, bottom, *_ = self.outs
        top_side_draw = self.top_side_draw
        bottom_side_draw = self.bottom_side_draw
        equations = []
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {}
        S = self.K * self.B
        for i in self.outs: eq_overall[i] = ones
        for i in process_inlets: eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        )
        
        # Top to bottom flows
        eq_outs = {}
        eq_outs[top] = ones
        eq_outs[bottom] = -S
        equations.append(
            (eq_outs, zeros)
        )
        
        # Top split flows
        if top_side_draw:
            eq_top_split = {
                top_side_draw: top_split * ones,
                top: top_split + minus_ones,
            }
            equations.append(
                (eq_top_split, zeros)
            )
        
        # Bottom split flows
        if bottom_side_draw:
            eq_bottom_split = {
                bottom_side_draw: bottom_split * ones,
                bottom: bottom_split + minus_ones,
            }
            equations.append(
                (eq_bottom_split, zeros)
            )
        
        return equations
    
    def _create_linear_equations(self, variable):
        # list[dict[Unit|Stream, float]]
        phases = self.phases
        partition = self.partition
        chemicals = self.chemicals
        pIDs = partition.IDs
        IDs = chemicals.IDs
        if pIDs != IDs:
            partition.IDs = IDs
            K = np.ones(chemicals.size)
            index = [IDs.index(i) for i in pIDs]
            for i, j in zip(index, partition.K): K[i] = j
            partition.K = K
            if phases == ('L', 'l') and not getattr(self, 'coupled_KL', None):
                if partition.gamma_y is not None:
                    gamma_y = np.ones(chemicals.size)
                    for i, j in zip(index, partition.gamma_y): gamma_y[i] = j
                    partition.gamma_y = gamma_y
        if (variable == 'K' and phases == ('g', 'l')
            or variable == 'K-pseudo' and phases == ('L', 'l')):
            if phases == ('g', 'l'):
                partition._run_decoupled_KTvle()
            elif phases == ('L', 'l'):
                partition._run_decoupled_Kgamma()
            else:
                raise NotImplementedError(f'K for phases {phases} is not yet implemented')
            eqs = [({self: np.ones(chemicals.size)}, partition.K)]
        elif variable == 'T':
            if phases == ('g', 'l'):
                eqs = []
            elif phases == ('L', 'l'):
                eqs = [self._create_temperature_departure_equation()]
            else:
                raise NotImplementedError(f'K for phases {phases} is not yet implemented')
        elif variable == 'B':
            if phases == ('g', 'l') and self.B_specification is None:
                eqs = [self._create_phase_ratio_departure_equation()]
            else:
                eqs = []
        elif variable == 'L' and phases == ('L', 'l'):
            if not getattr(self, 'coupled_KL', None):
                self.partition._run_decoupled_B()
            eqs = []
        elif variable == 'mol':
            eqs = self._create_mass_balance_equations()
        elif variable == 'mol-LLE' and phases == ('L', 'l'):
            if getattr(self, 'coupled_KL', None):
                eqs = []
            else:
                eqs = self._create_lle_mass_balance_equations()
        elif variable == 'KL' and getattr(self, 'coupled_KL', None):
            self.partition._run_lle(update=False)
            eqs = []
        else:
            eqs = []
        return eqs
    
    def _get_decoupled_variable(self, variable):
        if variable == 'T': 
            return self.T
        elif variable == 'L' and self.phases == ('L', 'l'): 
            if getattr(self, 'coupled_KL', None): return
            return self.B
        elif (variable == 'K-pseudo' and self.phases == ('L', 'l')
              or variable == 'K' and self.phases == ('g', 'l')):
            if getattr(self, 'coupled_KL', None): return
            return self.K
    
    def _update_decoupled_variable(self, variable, value):
        if variable == 'T':
            self.T = T = self.T + value
            for i in self.outs: i.T = T
        elif variable == 'B':
            self.B += value
        elif variable == 'L':
            self.B = value
        elif (variable == 'K-pseudo' or variable == 'K'): 
            self.K = value
    
# %%

class PhasePartition(Unit):
    _N_ins = 1
    _N_outs = 2
    strict_infeasibility_check = False
    
    def _init(self, phases, partition_data, top_chemical=None):
        self.partition_data = partition_data
        self.phases = phases
        self.top_chemical = top_chemical
        self.gamma_y = None
        self.IDs = None
        self.K = None
        self.B = None
        self.T = None
        self.Q = 0.
        self.B_specification = None
        self.B_fallback = 1
        for i, j in zip(self.outs, self.phases): i.phase = j 
        
    def _get_mixture(self, linked=True):
        if linked:
            try:
                ms = self._linked_multistream 
            except:
                outs = self.outs
                self._linked_multistream = ms = tmo.MultiStream.from_streams(outs)
            ms.copy_like(self.feed)
            return ms
        else:
            try:
                ms = self._unlinked_multistream
                ms.copy_like(self.feed)
            except:
                self._unlinked_multistream = ms = self.feed.copy()
                ms.phases = self.phases
            return ms
    
    def _get_arrays(self):
        if self.gamma_y is None:
            return {'K': self.K}
        else:
            return {'K': self.K, 'gamma_y': self.gamma_y}
    
    def _set_arrays(self, IDs, **kwargs):
        IDs_last = self.IDs
        IDs = tuple(IDs)
        if IDs_last and IDs_last != IDs and len(IDs_last) > len(IDs):
            size = len(IDs_last)
            index = [IDs_last.index(i) for i in IDs]
            for name, array in kwargs.items():
                last = getattr(self, name)
                if last.size != size:
                    last = np.ones(size)
                    setattr(self, name, last)
                for i, j in enumerate(index):
                    last[j] = array[i]
        else:
            for i, j in kwargs.items(): setattr(self, i, j)
            self.IDs = IDs
    
    def _get_activity_model(self):
        chemicals = self.chemicals
        index = chemicals.get_lle_indices(sum([i.mol for i in self.ins]).nonzero_keys())
        chemicals = chemicals.tuple
        lle_chemicals = [chemicals[i] for i in index]
        return self.thermo.Gamma(lle_chemicals), [i.ID for i in lle_chemicals], index
    
    def _run_decoupled_Kgamma(self, P=None): # Psuedo-equilibrium
        top, bottom = self.outs
        f_gamma, IDs, index = self._get_activity_model()
        T = self.T
        x = bottom.mol[index]
        x_sum = x.sum()
        if x_sum:
            x /= x_sum
        else:
            x = np.ones(x.size) / x.size
        gamma_x = f_gamma(x, T)
        gamma_y = self.gamma_y
        try:
            init_gamma = gamma_y is None or gamma_y.size != len(index)
        except:
            init_gamma = True
        if init_gamma:
            y = top.mol[index]
            y_sum = y.sum()
            if y_sum: 
                y /= y_sum
            else:
                y = np.ones(y.size) / y.size
            self.gamma_y = gamma_y = f_gamma(y, T)
        K = self.K
        K = gamma_x / gamma_y 
        y = K * x
        y /= y.sum()
        gamma_y = f_gamma(y, T)
        K = gamma_x / gamma_y
        good = (x != 0) | (y != 0)
        if not good.all():
            index, = np.where(good)
            IDs = [IDs[i] for i in index]
            gamma_y = [gamma_y[i] for i in index]
            K = [K[i] for i in index]
        self._set_arrays(IDs, gamma_y=gamma_y, K=K)
        
    def _run_decoupled_B(self, stacklevel=1): # Flash Rashford-Rice
        ms = self.feed.copy()
        ms.phases = self.phases
        top, bottom = ms
        data = self.partition_data
        try:
            if data and 'K' in data:
                phi = sep.partition(
                    ms, top, bottom, self.IDs, data['K'], 0.5, 
                    data.get('extract_chemicals') or data.get('top_chemicals'),
                    data.get('raffinate_chemicals') or data.get('bottom_chemicals'),
                    self.strict_infeasibility_check, stacklevel+1
                )
            else:
                phi = sep.partition(
                    ms, top, bottom, self.IDs, self.K, 0.5, 
                    None, None, self.strict_infeasibility_check,
                    stacklevel+1
                )
        except: 
            return self.B_fallback
        if phi <= 0 or phi >= 1: return
        self.B = phi / (1 - phi)
    
    def _run_decoupled_KTvle(self, P=None): # Bubble point
        top, bottom = self.outs
        if bottom.isempty():
            if top.isempty(): return
            p = top.dew_point_at_P(P)
        elif top.isempty():
            return
        else:
            p = bottom.bubble_point_at_P(P)
        # TODO: Note that solution decomposition method is bubble point
        x = p.x
        x[x == 0] = 1.
        K_new = p.y / p.x
        IDs = p.IDs
        self.T = p.T
        self._set_arrays(IDs, K=K_new)
    
    def _run_lle(self, P=None, update=True, top_chemical=None):
        if top_chemical is None: top_chemical = self.top_chemical
        else: self.top_chemical = top_chemical
        ms = self._get_mixture(update)
        eq = ms.lle
        data = self.partition_data
        if data and 'K' in data:
            ms.phases = self.phases
            top, bottom = ms
            phi = sep.partition(
                ms, top, bottom, self.IDs, data['K'], 0.5, 
                data.get('extract_chemicals') or data.get('top_chemicals'),
                data.get('raffinate_chemicals') or data.get('bottom_chemicals'),
                self.strict_infeasibility_check,1
            )
            if phi == 1:
                self.B = np.inf
            else:
                self.B = phi / (1 - phi)
        else:
            if update:
                eq(T=ms.T, P=P, top_chemical=top_chemical, update=update)
                lle_chemicals, K_new, phi = eq._lle_chemicals, eq._K, eq._phi
            else:
                lle_chemicals, K_new, phi = eq(T=ms.T, P=P, top_chemical=top_chemical, update=update)
            if phi == 1 or phi is None:
                self.B = np.inf
                return
            else:
                self.B = phi / (1 - phi)
            self.T = ms.T
            IDs = tuple([i.ID for i in lle_chemicals])
            self._set_arrays(IDs, K=K_new)
    
    def _run_vle(self, P=None, update=True):
        ms = self._get_mixture(update)
        B = self.B_specification
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
        ms.vle(P=P or ms.P, H=H, V=V)
        index = ms.vle._index
        IDs = ms.chemicals.IDs
        IDs = tuple([IDs[i] for i in index])
        L_mol = ms.imol['l', IDs]
        L_total = L_mol.sum()
        if L_total: 
            x_mol = L_mol / L_total
            x_mol[x_mol == 0] = 1e-9
        else:
            x_mol = 1
        V_mol = ms.imol['g', IDs]
        V_total = V_mol.sum()
        if V_total: 
            y_mol = V_mol / V_total
        else:
            y_mol = 0
        K_new = y_mol / x_mol
        if B is None: 
            if not L_total:
                self.B = inf
            else:
                self.B = V_total / L_total
        self.T = ms.T
        self._set_arrays(IDs, K=K_new)
    
    def _run(self):
        if self.phases == ('g', 'l'):
            self._run_vle()
        else:
            self._run_lle()


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
    top_chemical : str
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
    >>> extract, raffinate = MSE.outs
    >>> extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.83
    >>> extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    0.99
    >>> raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    0.82
    
    Simulate 10-stage extraction with user defined partition coefficients:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> import numpy as np
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=10, ins=[feed, solvent], phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.451e-01, 1.380e+00, 2.958e+03]),
    ...         'IDs': ('Water', 'Methanol', 'Octanol'),
    ...         'phi': 0.5899728891780545, # Initial phase fraction guess. This is optional.
    ...     }
    ... )
    >>> extract, raffinate = MSE.outs
    >>> MSE.simulate()
    >>> extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.99
    >>> extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    0.99
    >>> raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    0.82
    
    Because octanol and water do not mix well, it may be a good idea to assume
    that these solvents do not mix at all:
        
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> import numpy as np
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=20, ins=[feed, solvent], phases=('L', 'l'),
    ...     partition_data={
    ...         'K': np.array([1.38]),
    ...         'IDs': ('Methanol',),
    ...         'raffinate_chemicals': ('Water',),
    ...         'extract_chemicals': ('Octanol',),
    ...     }
    ... )
    >>> MSE.simulate()
    >>> extract, raffinate = MSE.outs
    >>> extract.imol['Methanol'] / feed.imol['Methanol'] # Recovery
    0.99
    >>> extract.imol['Octanol'] / solvent.imol['Octanol'] # Solvent stays in extract
    1.0
    >>> raffinate.imol['Water'] / feed.imol['Water'] # Carrier remains in raffinate
    1.0
       
    Simulate with a feed at the 4th stage:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> import numpy as np
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
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
    >>> extract, raffinate = MSE.outs
    >>> extract.imol['Methanol'] / (feed.imol['Methanol'] + dilute_feed.imol['Methanol']) # Recovery
    0.93
    
    Simulate with a 60% extract side draw at the 2nd stage:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Methanol', 'Octanol'], cache=True)
    >>> import numpy as np
    >>> feed = bst.Stream('feed', Water=5000, Methanol=500)
    >>> solvent = bst.Stream('solvent', Octanol=5000)
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
    0.92
    
    Simulate stripping column with 2 stages
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['AceticAcid', 'EthylAcetate', 'Water', 'MTBE'], cache=True)
    >>> feed = bst.Stream('feed', Water=75, AceticAcid=5, MTBE=20, T=320)
    >>> steam = bst.Stream('steam', Water=100, phase='g', T=390)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=2, ins=[feed, steam], feed_stages=[0, -1],
    ...     outs=['vapor', 'liquid'],
    ...     phases=('g', 'l'),
    ... )
    >>> MSE.simulate()
    >>> vapor, liquid = MSE.outs
    >>> vapor.imol['MTBE'] / feed.imol['MTBE']
    0.99
    >>> vapor.imol['Water'] / (feed.imol['Water'] + steam.imol['Water'])
    0.42
    >>> vapor.imol['AceticAcid'] / feed.imol['AceticAcid']
    0.74
    
    Simulate distillation column with 5 stages, a 0.673 reflux ratio, 
    2.57 boilup ratio, and feed at stage 2:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Ethanol=80, Water=100, T=80.215 + 273.15)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed], feed_stages=[2],
    ...     outs=['vapor', 'liquid'],
    ...     stage_specifications={0: ('Reflux', 0.673), -1: ('Boilup', 2.57)},
    ...     phases=('g', 'l'),
    ... )
    >>> MSE.simulate()
    >>> vapor, liquid = MSE.outs
    >>> vapor.imol['Ethanol'] / feed.imol['Ethanol']
    0.96
    >>> vapor.imol['Ethanol'] / vapor.F_mol
    0.69
    
    Simulate the same distillation column with a full condenser, 5 stages, a 0.673 reflux ratio, 
    2.57 boilup ratio, and feed at stage 2:
    
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Ethanol=80, Water=100, T=80.215 + 273.15)
    >>> MSE = bst.MultiStageEquilibrium(N_stages=5, ins=[feed], feed_stages=[2],
    ...     outs=['vapor', 'liquid', 'distillate'],
    ...     stage_specifications={0: ('Reflux', float('inf')), -1: ('Boilup', 2.57)},
    ...     bottom_side_draws={0: 0.673 / (1 + 0.673)}
    ... )
    >>> MSE.simulate()
    >>> vapor, liquid, distillate = MSE.outs
    >>> distillate.imol['Ethanol'] / feed.imol['Ethanol']
    0.81
    >>> distillate.imol['Ethanol'] / distillate.F_mol
    0.70
    
    """
    _N_ins = 2
    _N_outs = 2
    default_maxiter = 5
    default_max_attempts = 10
    default_fallback_maxiter = 1
    default_molar_tolerance = 0.1
    default_relative_molar_tolerance = 0.001
    default_algorithm = 'root'
    available_algorithms = {'root', 'optimize'}
    default_methods = {
        'root': 'fixed-point',
        'optimize': 'SLSQP',
        'SurPASS': 'differential evolution', 
    }
    
    #: Method definitions for convergence
    root_options: dict[str, tuple[Callable, bool, dict]] = {
        'fixed-point': (flx.conditional_fixed_point, True, {}),
    }
    optimize_options: dict[str, tuple[Callable, dict]] = {
        'SLSQP': (minimize, {'tol': 1e-3, 'method': 'SLSQP'})
    }
    SurPASS_options: dict[str, tuple[Callable, dict]] = {
        'differential evolution': (
            differential_evolution, {'seed': 0, 'popsize': 12, 'tol': 1e-6}
        ),
    }
    auxiliary_unit_names = (
        'stages',
    )
    _side_draw_names = ('top_side_draws', 'bottom_side_draws')
    
    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(cls, *args, **kwargs)
        if '_side_draw_names' in cls.__dict__:
            top, bottom = cls._side_draw_names
            setattr(
                cls, top, 
                property(
                    lambda self: self.top_side_draws,
                    lambda self, value: setattr(self, 'top_side_draws', value)
                )
            )
            setattr(
                cls, bottom, 
                property(
                    lambda self: self.bottom_side_draws,
                    lambda self, value: setattr(self, 'bottom_side_draws', value)
                )
            )
    
    def __init__(self,  ID='', ins=None, outs=(), thermo=None, **kwargs):
        if 'feed_stages' in kwargs: self._N_ins = len(kwargs['feed_stages'])
        top_side_draws, bottom_side_draws = self._side_draw_names
        N_outs = 2
        if top_side_draws in kwargs: N_outs += len(kwargs[top_side_draws]) 
        if bottom_side_draws in kwargs: N_outs += len(kwargs[bottom_side_draws]) 
        self._N_outs = N_outs
        Unit.__init__(self, ID, ins, outs, thermo, **kwargs)
    
    def _init(self,
            N_stages, 
            top_side_draws=None,
            bottom_side_draws=None, 
            feed_stages=None, 
            phases=None, 
            P=101325, 
            stage_specifications=None, 
            partition_data=None, 
            top_chemical=None, 
            use_cache=None,
            collapsed_init=True,
            algorithm=None,
            method=None,
            maxiter=None,
            inside_out=None,
        ):
        # For VLE look for best published algorithm (don't try simple methods that fail often)
        if phases is None: phases = ('g', 'l')
        if feed_stages is None: feed_stages = (0, -1)
        if stage_specifications is None: stage_specifications = {}
        elif not isinstance(stage_specifications, dict): stage_specifications = dict(stage_specifications)
        if top_side_draws is None: top_side_draws = {}
        elif not isinstance(top_side_draws, dict): top_side_draws = dict(top_side_draws)
        if bottom_side_draws is None: bottom_side_draws = {}
        elif not isinstance(bottom_side_draws, dict): bottom_side_draws = dict(bottom_side_draws)
        if partition_data is None: partition_data = {}
        self.multi_stream = tmo.MultiStream(None, P=P, phases=phases, thermo=self.thermo)
        self.N_stages = N_stages
        self.P = P
        self.phases = phases = self.multi_stream.phases # Corrected order
        self._has_vle = 'g' in phases
        self._has_lle = 'L' in phases
        top_mark = 2 + len(top_side_draws)
        tsd_iter = iter(self.outs[2:top_mark])
        bsd_iter = iter(self.outs[top_mark:])
        last_stage = None
        self._top_split = top_splits = np.zeros(N_stages)
        self._bottom_split = bottom_splits = np.zeros(N_stages)
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
                outs.append(None)
            if i == N_stages - 1: 
                outs.append(
                    self-1 # raffinate or liquid
                )
            else:
                outs.append(None)
            if i in top_side_draws:
                outs.append(next(tsd_iter))
                top_split = top_side_draws[i]
                top_splits[i] = top_split 
            else: 
                top_split = 0
            if i in bottom_side_draws:
                outs.append(next(bsd_iter))
                bottom_split = bottom_side_draws[i]
                bottom_splits[i] = bottom_split
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
        self._asplit_left = 1 - top_splits
        self._bsplit_left = 1 - bottom_splits
        self._asplit_1 = top_splits - 1
        self._bsplit_1 = bottom_splits - 1
        self.partitions = [i.partition for i in stages]
        self.top_chemical = top_chemical
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
        self.maxiter = self.default_maxiter if maxiter is None else maxiter
        
        #: [int] Maximum number of iterations for fallback algorithm.
        self.fallback_maxiter = self.default_fallback_maxiter

        #: [float] Molar tolerance (kmol/hr)
        self.molar_tolerance = self.default_molar_tolerance

        #: [float] Relative molar tolerance
        self.relative_molar_tolerance = self.default_relative_molar_tolerance
        
        self.use_cache = True if use_cache else False
        
        self.collapsed_init = collapsed_init
        
        self.algorithm = self.default_algorithm if algorithm is None else algorithm
        
        self.method = self.default_methods[self.algorithm] if method is None else method
    
        self.inside_out = inside_out
        
        self.max_attempts = self.default_max_attempts
    
    def correct_overall_mass_balance(self):
        outmol = sum([i.mol for i in self.outs])
        inmol = sum([i.mol for i in self.ins])
        try:
            factor = inmol / outmol
        except:
            pass
        else:
            for i in self.outs: i.mol *= factor
    
    def material_errors(self):
        errors = []
        stages = self.stages
        IDs = self.multi_stream.chemicals.IDs
        for stage in stages:
            errors.append(
                sum([i.imol[IDs] for i in stage.ins],
                    -sum([i.imol[IDs] for i in stage.outs]))
            )
        return pd.DataFrame(errors, columns=IDs)
    
    def set_flow_rates(self, top_flows):
        top, bottom = self.multi_stream.phases
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        index = self._update_index
        top_flows[top_flows < 0] = 0
        bottom_flows = mass_balance(
            top_flows, self.feed_flows, self._asplit_left, self._bsplit_left, 
            np.zeros(N_stages, bool), self.N_stages, self._N_chemicals
        )
        bottom_flows[bottom_flows < 0] = 0
        for i in range_stages:
            stage = stages[i]
            partition = stage.partition
            s_top, s_bot = partition.outs
            s_top.mol[index] = top_flows[i]
            s_bot.mol[index] = bottom_flows[i]
            for i in stage.splitters: i._run()
        
    def set_flow_rates_old(self, top_flows):
        top, bottom = self.multi_stream.phases
        flow_tol = -1e-6 * self.multi_stream.mol
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        index = self._update_index
        top_flows[top_flows < 0.] = 0.
        has_infeasible_flow = True
        infeasible_checks = set()
        while has_infeasible_flow:
            has_infeasible_flow = False
            for i in range_stages:
                stage = stages[i]
                partition = stage.partition
                s_top, _ = partition.outs
                s_top.mol[index] = top_flows[i]
                if stage.top_split: stage.splitters[0]._run()
            for i in range_stages:
                stage = stages[i]
                partition = stage.partition
                s_top, s_bottom = partition.outs
                bottom_flow = sum([i.mol for i in stage.ins]) - s_top.mol
                mask = bottom_flow < 0.
                if mask.any():
                    has_infeasible_flow = (bottom_flow[mask] < flow_tol[mask]).any()
                    if i not in infeasible_checks and has_infeasible_flow:
                        infeasible_checks.add(i)
                        infeasible_index, = np.where(mask[index])
                        # TODO: Find algebraic solution to keeping top flow rates within feasible region.
                        # This is only a temporary solution.
                        infeasible_flow = bottom_flow[mask]
                        top_flows[i, infeasible_index] += infeasible_flow
                        break
                    else:
                        has_infeasible_flow = False
                        bottom_flow[mask] = 0.
                s_bottom.mol[:] = bottom_flow
                if stage.bottom_split: stage.splitters[-1]._run()
        # self.correct_overall_mass_balance()
            
    def _run(self):
        if all([i.isempty() for i in self.ins]): 
            for i in self.outs: i.empty()
            return
        top_flow_rates = self.hot_start()
        algorithm = self.algorithm
        if algorithm == 'root':
            self.converged = True
            for i in range(self.max_attempts):
                self.attempt = i
                self.iter = 0
                self.fallback_iter = 0
                method = self.method
                solver, conditional, options = self.root_options[method]
                if conditional:
                    top_flow_rates = solver(self._conditional_iter, top_flow_rates)
                else:
                    try:
                        result = solver(self._root_iter, self.get_KTBs().flatten(), **options)
                        break
                    except RuntimeError: # Fall back to fixed-point
                        method = 'fixed-point'
                        solver, conditional, options = self.root_options[method]
                        top_flow_rates = solver(self._conditional_iter, top_flow_rates)
                if method == 'fixed-point' and self.iter == self.maxiter:
                    top_flow_rates = flx.conditional_fixed_point(
                        self._sequential_iter, 
                        top_flow_rates,
                    )
                    if self.fallback_iter < self.fallback_maxiter:
                        break
                else:
                    break
            else:
                self.converged = False
        elif algorithm == 'optimize':
            solver, options = self.optimize_options[self.method]
            self.constraints = constraints = []
            stages = self.stages
            m, n = self.N_stages, self._N_chemicals
            last_stage = m - 1
            feed_flows, asplit_1, bsplit_1, _ = self._iter_args
            for i, stage in enumerate(stages):
                if i == 0:
                    args = (i,)
                    f = lambda x, i: feed_flows[i] - x[(i+1)*n:(i+2)*n] * asplit_1[i+1] - x[i*n:(i+1)*n] + 1e-6
                elif i == last_stage:
                    args_last = args
                    args = (i, f, args_last)
                    f = lambda x, i, f, args_last: feed_flows[i] + f(x, *args_last) - x[i*n:] + 1e-6
                else:
                    args_last = args
                    args = (i, f, args_last)
                    f = lambda x, i, f, args_last: feed_flows[i] + f(x, *args_last) - x[(i+1)*n:(i+2)*n] * asplit_1[i+1] - x[i*n:(i+1)*n] + 1e-6
                constraints.append(
                    dict(type='ineq', fun=f, args=args)
                )
            result = minimize(
                self._overall_error, 
                self.get_top_flow_rates_flat(),
                constraints=constraints,
                bounds=[(0, None)] * (m * n),
                **options,
            )
            self.set_flow_rates(result.x.reshape([m, n]))
        else:
            raise RuntimeError(
                f'invalid algorithm {algorithm!r}, only {self.available_algorithms} are allowed'
            )
        self.correct_overall_mass_balance()
    
    def _hot_start_phase_ratios_iter(self, 
            top_flow_rates, *args
        ):
        bottom_flow_rates = hot_start_bottom_flow_rates(
            top_flow_rates, *args
        )
        top_flow_rates = hot_start_top_flow_rates(
            bottom_flow_rates, *args
        )
        return top_flow_rates
        
    def hot_start_phase_ratios(self):
        stages = self.stages
        stage_index = []
        phase_ratios = []
        for i in list(self.stage_specifications):
            B = stages[i].partition.B_specification
            if B is None: continue 
            stage_index.append(i)
            phase_ratios.append(B)
        stage_index = np.array(stage_index, dtype=int)
        phase_ratios = np.array(phase_ratios, dtype=float)
        feeds = self.ins
        feed_stages = self.feed_stages
        top_feed_flows = 0 * self.feed_flows
        bottom_feed_flows = top_feed_flows.copy()
        top_flow_rates = top_feed_flows.copy()
        index = self._update_index
        for feed, stage in zip(feeds, feed_stages):
            if len(feed.phases) > 1 and 'g' in feed.phases:
                top_feed_flows[stage, :] += feed['g'].mol[index]
            elif feed.phase != 'g':
                continue
            else:
                top_feed_flows[stage, :] += feed.mol[index]
        for feed, stage in zip(feeds, feed_stages):
            if len(feed.phases) > 1 and 'g' not in feed.phases:
                bottom_feed_flows[stage, :] += feed['l'].mol[index]
            elif feed.phase == 'g': 
                continue
            else:
                bottom_feed_flows[stage, :] += feed.mol[index]
        feed_flows, asplit_1, bsplit_1, N_stages = self._iter_args
        args = (
            phase_ratios, np.array(stage_index), top_feed_flows,
            bottom_feed_flows, asplit_1, bsplit_1, N_stages
        )
        top_flow_rates = flx.wegstein(
            self._hot_start_phase_ratios_iter,
            top_flow_rates, args=args, xtol=self.relative_molar_tolerance,
            checkiter=False,
        )
        bottom_flow_rates = hot_start_bottom_flow_rates(
            top_flow_rates, *args
        )
        bf = bottom_flow_rates.sum(axis=1)
        bf[bf == 0] = 1e-32
        return top_flow_rates.sum(axis=1) / bf
    
    def hot_start_collapsed_stages(self,
            all_stages, feed_stages, stage_specifications,
            top_side_draws, bottom_side_draws,
        ):
        N_stages = len(all_stages)
        stage_map = {j: i for i, j in enumerate(sorted(all_stages))}
        feed_stages = [stage_map[i] for i in feed_stages]
        stage_specifications = {stage_map[i]: j for i, j in stage_specifications.items()}
        top_side_draws = {stage_map[i]: j for i, j in top_side_draws.items()}
        bottom_side_draws = {stage_map[i]: j for i, j in bottom_side_draws.items()}
        self.collapsed = collapsed = MultiStageEquilibrium(
            '.collapsed', 
            ins=[i.copy() for i in self.ins],
            outs=[i.copy() for i in self.outs],
            N_stages=N_stages,
            feed_stages=feed_stages,
            stage_specifications=stage_specifications,
            phases=self.multi_stream.phases,
            top_side_draws=top_side_draws,
            bottom_side_draws=bottom_side_draws,  
            P=self.P, 
            partition_data=self.partition_data,
            top_chemical=self.top_chemical, 
            use_cache=self.use_cache,
            thermo=self.thermo
        )
        collapsed._run()
        collapsed_stages = collapsed.stages
        partitions = self.partitions
        stages = self.stages
        for i in range(self.N_stages):
            if i in all_stages:
                collapsed_partition = collapsed_stages[stage_map[i]].partition
                partition = partitions[i]
                partition.T = collapsed_partition.T
                partition.B = collapsed_partition.B
                T = collapsed_partition.T
                for i in partition.outs + stages[i].outs: i.T = T 
                partition.K = collapsed_partition.K
                partition.gamma_y = collapsed_partition.gamma_y
        self.interpolate_missing_variables()
                
    def hot_start(self):
        ms = self.multi_stream
        feeds = self.ins
        feed_stages = self.feed_stages
        stages = self.stages
        partitions = self.partitions
        N_stages = self.N_stages
        chemicals = self.chemicals
        top_phase, bottom_phase = ms.phases
        eq = 'vle' if top_phase == 'g' else 'lle'
        ms.mix_from(feeds)
        ms.P = self.P
        if eq == 'lle':
            self.top_chemical = top_chemical = self.top_chemical or feeds[1].main_chemical
            for i in partitions: i.top_chemical = top_chemical
        data = self.partition_data
        if data:
            top_chemicals = data.get('extract_chemicals') or data.get('vapor_chemicals', [])
            bottom_chemicals = data.get('raffinate_chemicals') or data.get('liquid_chemicals', [])
            for i in chemicals.light_chemicals:
                i = i.ID
                if i in top_chemicals or i in bottom_chemicals: continue
                top_chemicals.append(i)
            for i in chemicals.heavy_chemicals:
                i = i.ID
                if i in top_chemicals or i in bottom_chemicals: continue
                bottom_chemicals.append(i)
        else:
            top_chemicals = [i.ID for i in chemicals.light_chemicals]
            bottom_chemicals = [i.ID for i in chemicals.heavy_chemicals]
        if eq == 'lle':
            IDs = data['IDs'] if 'IDs' in data else [i.ID for i in ms.lle_chemicals]
        else:
            IDs = data['IDs'] if 'IDs' in data else [i.ID for i in ms.vle_chemicals]
        IDs = tuple(IDs)
        self._N_chemicals = N_chemicals = len(IDs)
        self._update_index = index = ms.chemicals.get_index(IDs)
        self.feed_flows = feed_flows = np.zeros([N_stages, N_chemicals])
        self.feed_enthalpies = feed_enthalpies = np.zeros(N_stages)
        for feed, stage in zip(feeds, feed_stages):
            feed_flows[stage, :] += feed.mol[index]
            feed_enthalpies[stage] = feed.H
        self._iter_args = (feed_flows, self._asplit_1, self._bsplit_1, self.N_stages)
        feed_stages = [(i if i >= 0 else N_stages + i) for i in self.feed_stages]
        stage_specifications = {(i if i >= 0 else N_stages + i): j for i, j in self.stage_specifications.items()}
        top_side_draws = {(i if i >= 0 else N_stages + i): j for i, j in self.top_side_draws.items()}
        bottom_side_draws = {(i if i >= 0 else N_stages + i): j for i, j in self.bottom_side_draws.items()}
        self.key_stages = key_stages = set([*feed_stages, *stage_specifications, *top_side_draws, *bottom_side_draws])
        if (self.use_cache 
            and all([i.IDs == IDs for i in partitions])): # Use last set of data
            pass
        elif self.collapsed_init and len(key_stages) != self.N_stages:
            self.hot_start_collapsed_stages(
                key_stages, feed_stages, stage_specifications,
                top_side_draws, bottom_side_draws,
            )
        else:
            if data and 'K' in data: 
                top, bottom = ms
                K = data['K']
                phi = data.get('phi') or top.imol[IDs].sum() / ms.imol[IDs].sum()
                data['phi'] = phi = sep.partition(ms, top, bottom, IDs, K, phi,
                                                  top_chemicals, bottom_chemicals)
                B = inf if phi == 1 else phi / (1 - phi)
                T = ms.T
                for i in partitions: 
                    if i.B_specification is None: i.B = B
                    i.T = T
                    
            elif eq == 'lle':
                lle = ms.lle
                T = ms.T
                lle(T, top_chemical=top_chemical)
                K = lle._K
                phi = lle._phi
                B = inf if phi == 1 else phi / (1 - phi)
                y = ms.imol['L', IDs]
                y /= y.sum()
                f_gamma = self.thermo.Gamma([chemicals[i] for i in IDs])
                gamma_y = f_gamma(y, T)
                for i in partitions: 
                    i.B = B
                    i.T = T
                    i.gamma_y = gamma_y
                    for j in i.outs: j.T = T
            else:
                P = self.P
                if self.stage_specifications:
                    dp = ms.dew_point_at_P(P=P, IDs=IDs)
                    T_bot = dp.T
                    bp = ms.bubble_point_at_P(P=P, IDs=IDs)
                    T_top = bp.T
                    dT_stage = (T_bot - T_top) / N_stages
                    phase_ratios = self.hot_start_phase_ratios()
                    K = bp.y / bp.z
                    for i, B in enumerate(phase_ratios):
                        partition = partitions[i]
                        if partition.B_specification is None: partition.B = B
                        partition.T = T = T_top + i * dT_stage
                        for s in partition.outs: s.T = T
                else:
                    vle = ms.vle
                    vle(H=ms.H, P=P)
                    L_mol = ms.imol['l', IDs]
                    L_mol_net = L_mol.sum()
                    if L_mol_net: x_mol = L_mol / L_mol.sum()
                    else: x_mol = np.ones(N_chemicals, float) / N_chemicals
                    V_mol = ms.imol['g', IDs]
                    y_mol = V_mol / V_mol.sum()
                    K = y_mol / x_mol
                    phi = ms.V
                    B = 1 / (1 - phi)
                    T = ms.T
                    for partition in partitions:
                        partition.T = T
                        partition.B = B
                        for i in partition.outs: i.T = T
            for i in partitions: i.K = K
            N_chemicals = len(index)
        if top_chemicals:
            top_side_draws = self.top_side_draws
            n = len(top_chemicals)
            b = np.ones([N_stages, n])
            c = self._asplit_1[1:]
            d = np.zeros([N_stages, n])
            for feed, stage in zip(feeds, feed_stages):
                d[stage] += feed.imol[top_chemicals]
            top_flow_rates = solve_RBDMA(b, c, d)
            for partition, flows in zip(partitions, top_flow_rates):
                partition.outs[0].imol[top_chemicals] = flows
        if bottom_chemicals:
            bottom_side_draws = self.bottom_side_draws
            a = self._bsplit_1[:-1]
            n = len(bottom_chemicals)
            b = np.ones([N_stages, n])
            d = np.zeros([N_stages, n])
            for feed, stage in zip(feeds, feed_stages):
                d[stage] += feed.imol[bottom_chemicals]
            bottom_flow_rates = solve_LBDMA(a, b, d)
            for partition, b in zip(partitions, bottom_flow_rates):
                partition.outs[1].imol[bottom_chemicals] = b
        if top_chemicals or bottom_chemicals:
            for i in stages:
                for s in i.splitters: s._run()
        for i in partitions: i.IDs = IDs
        return self.run_mass_balance()
    
    # TODO: This method is not working well. 
    # It should be mathematically the same as the departure method,
    # but it gives results that are extremely off.
    # def get_energy_balance_temperatures(self):
    #     # ENERGY BALANCE
    #     # Hv1 + Cv1*dT1 + Hl1 + Cl1*dT1 - Hv2 - Cv2*dT2 - Hl0 - Cl0*T0 = Q1
    #     # (Cv1 + Cl1)dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - Hv1 - Hl1 + Hv2 + Hl0
    #     # C1*T1 - Cv2*T2 - Cl0*T0 = Q1 - H_out + H_in + C1*T1ref - Cv2*T2ref - Cl0*T0ref
    #     stages = self.stages
    #     N_stages = self.N_stages
    #     a = np.zeros(N_stages)
    #     b = a.copy()
    #     c = a.copy()
    #     d = a.copy()
    #     stage_mid = stages[0]
    #     stage_bottom = stages[1]
    #     partition_mid = stage_mid.partition
    #     partition_bottom = stage_bottom.partition
    #     C_out = sum([i.C for i in partition_mid.outs])
    #     C_bottom = stage_bottom.outs[0].C
    #     b[0] = C_out
    #     c[0] = -C_bottom
    #     Q = (partition_mid.Q or 0.) + C_out * partition_mid.T - C_bottom * partition_bottom.T
    #     if partition_mid.B_specification is not None:
    #         Q += stage_mid.H_in - partition_mid.H_out
    #     d[0] = Q
    #     for i in range(1, N_stages-1):
    #         stage_top = stage_mid
    #         stage_mid = stage_bottom
    #         stage_bottom = stages[i+1]
    #         partition_mid = partition_bottom
    #         partition_bottom = stage_bottom.partition
    #         C_out = sum([i.C for i in partition_mid.outs])
    #         C_top = stage_top.outs[1].C
    #         C_bottom = stage_bottom.outs[0].C
    #         a[i] = -C_top
    #         b[i] = C_out
    #         c[i] = -C_bottom
    #         Q = (partition_mid.Q or 0.) + C_out * partition_mid.T - C_bottom * partition_bottom.T - C_top * stage_top.partition.T
    #         if partition_mid.B_specification is not None:
    #             Q += stage_mid.H_in - partition_mid.H_out
    #         d[i] = Q
    #     stage_top = stage_mid
    #     stage_mid = stage_bottom
    #     partition_mid = partition_bottom
    #     C_out = sum([i.C for i in partition_mid.outs])
    #     C_top = stage_top.outs[1].C
    #     a[-1] = -C_top
    #     b[-1] = C_out
    #     Q = (partition_mid.Q or 0.) + C_out * partition_mid.T - C_top * stage_top.partition.T
    #     if partition_mid.B_specification is not None:
    #         Q += stage_mid.H_in - partition_mid.H_out
    #     d[-1] = Q
    #     return solve_TDMA(a, b, c, d)
    
    # Old slow algorithm, here for legacy purposes
    # def get_energy_balance_phase_ratio_departures(self):
    #     # ENERGY BALANCE
    #     # hv1*V1 + hl1*L1 - hv2*V2 - hl0*L0 = Q1
    #     # hV1*L1*B1 + hl1*L1 - hv2*L2*B2 - hl0*L0 = Q1
    #     # hV1*L1*B1 - hv2*L2*B2 = Q1 - Hl1 + Hl0
    #     # Hv1 + hV1*L1*dB1 - Hv2 - hv2*L2*dB2 = Q1 - Hl1 + Hl0
    #     # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
    #     stages = self.stages
    #     N_stages = self.N_stages
    #     b = np.zeros(N_stages)
    #     c = b.copy()
    #     d = c.copy()
    #     for i in range(N_stages-1):
    #         stage_mid = stages[i]
    #         stage_bottom = stages[i+1]
    #         partition_mid = stage_mid.partition
    #         partition_bottom = stage_bottom.partition
    #         mid = partition_mid.B_specification is None
    #         bottom = partition_bottom.B_specification is None
    #         vapor, liquid = partition_mid.outs
    #         vapor_bottom, liquid_bottom = partition_bottom.outs
    #         if mid:
    #             if vapor.isempty():
    #                 liquid.phase = 'g'
    #                 b[i] = liquid.H
    #                 liquid.phase = 'l'
    #             else:
    #                 b[i] = vapor.h * liquid.F_mol
    #         if bottom:
    #             if vapor_bottom.isempty():
    #                 liquid_bottom.phase = 'g'
    #                 c[i] = liquid_bottom.H
    #                 liquid_bottom.phase = 'l'
    #             else:
    #                 c[i] = -vapor_bottom.h * liquid_bottom.F_mol
    #         Q = partition_mid.Q or 0.
    #         if mid: d[i] = Q + stage_mid.H_in - partition_mid.H_out
    #     if bottom:
    #         if vapor_bottom.isempty():
    #             liquid_bottom.phase = 'g'
    #             b[i] = liquid_bottom.H
    #             liquid_bottom.phase = 'l'
    #         else:
    #             b[-1] = vapor_bottom.h * liquid_bottom.F_mol
    #         Q = partition_bottom.Q or 0.
    #         d[-1] = Q + stage_bottom.H_in - partition_bottom.H_out
    #     return solve_RBDMA_1D_careful(b, c, d)
    
    # Old method, here for legacy purposes
    # def get_energy_balance_temperature_departures_old(self):
    #     # ENERGY BALANCE
    #     # Hv1 + Cv1*(dT1) + Hl1 + Cl1*dT1 - Hv2 - Cv2*dT2 - Hl0 - Cl0 = Q1
    #     # (Cv1 + Cl1)dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - Hv1 - Hl1 + Hv2 + Hl0
    #     # C1dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - H_out + H_in
    #     stages = self.stages
    #     N_stages = self.N_stages
    #     a = np.zeros(N_stages)
    #     b = a.copy()
    #     c = a.copy()
    #     d = a.copy()
    #     stage_mid = stages[0]
    #     stage_bottom = stages[1]
    #     partition_mid = stage_mid.partition
    #     partition_bottom = stage_bottom.partition
    #     C_out = sum([i.C for i in partition_mid.outs])
    #     C_bottom = stage_bottom.outs[0].C
    #     b[0] = C_out
    #     c[0] = -C_bottom
    #     Q = partition_mid.Q or 0.
    #     if partition_mid.B_specification is None: d[0] = Q + stage_mid.H_in - partition_mid.H_out
    #     for i in range(1, N_stages-1):
    #         stage_top = stage_mid
    #         stage_mid = stage_bottom
    #         stage_bottom = stages[i+1]
    #         partition_mid = partition_bottom
    #         partition_bottom = stage_bottom.partition
    #         C_out = sum([i.C for i in partition_mid.outs])
    #         C_top = stage_top.outs[1].C
    #         C_bottom = stage_bottom.outs[0].C
    #         a[i] = -C_top
    #         b[i] = C_out
    #         c[i] = -C_bottom
    #         Q = partition_mid.Q or 0.
    #         if partition_mid.B_specification is None: d[i] = Q + stage_mid.H_in - partition_mid.H_out
    #     stage_top = stage_mid
    #     stage_mid = stage_bottom
    #     partition_mid = partition_bottom
    #     C_out = sum([i.C for i in partition_mid.outs])
    #     C_top = stage_top.outs[1].C
    #     a[-1] = -C_top
    #     b[-1] = C_out
    #     Q = partition_mid.Q or 0.
    #     if partition_mid.B_specification is None: d[-1] = Q + stage_mid.H_in - partition_mid.H_out
    #     return solve_TDMA(a, b, c, d)
    
    def get_energy_balance_temperature_departures(self):
        partitions = self.partitions
        N_stages = self.N_stages
        Cl = np.zeros(N_stages)
        Cv = Cl.copy()
        Hv = Cl.copy()
        Hl = Cl.copy()
        for i, j in enumerate(partitions):
            top, bottom = j.outs
            Hl[i] = bottom.H
            Hv[i] = top.H
            Cl[i] = bottom.C
            Cv[i] = top.C
        try:
            return temperature_departures(
                Cv, Cl, Hv, Hl, self._asplit_left, self._bsplit_left,
                N_stages, self.feed_enthalpies
            )
        except:
            return np.zeros(N_stages)
    
    def get_energy_balance_phase_ratio_departures(self):
        # ENERGY BALANCE
        # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
        partitions = self.partitions
        N_stages = self.N_stages
        L = np.zeros(N_stages)
        V = L.copy()
        hv = L.copy()
        hl = L.copy()
        specification_index = []
        missing = []
        for i, j in enumerate(partitions):
            top, bottom = j.outs
            Li = bottom.F_mol
            Vi = top.F_mol
            L[i] = Li
            V[i] = Vi
            if Vi == 0:
                if Li == 0:  
                    hv[i] = None
                    hl[i] = None
                    if j.B_specification: specification_index.append(i)
                    missing.append(i)
                    continue
                bottom.phase = 'g'
                hv[i] = bottom.h
                bottom.phase = 'l'
            else:
                hv[i] = top.h
            if Li == 0:
                top.phase = 'l'
                hl[i] = bottom.h
                top.phase = 'g'
            else:
                hl[i] = bottom.h
            if j.B_specification: specification_index.append(i)
        if missing:
            neighbors = get_neighbors(missing=missing, size=N_stages)
            hv = fillmissing(neighbors, hv)
            hl = fillmissing(neighbors, hl)
        return phase_ratio_departures(
            L, V, hl, hv, 
            self._asplit_1, 
            self._asplit_left,
            self._bsplit_left,
            N_stages,
            np.array(specification_index, dtype=int),
            self.feed_enthalpies,
        )
        
    def update_energy_balance_phase_ratios(self):
        dBs = self.get_energy_balance_phase_ratio_departures()
        for i, dB in zip(self.partitions, dBs):
            if i.B_specification is None: i.B += dB
    
    def update_energy_balance_temperatures(self):
        dTs = self.get_energy_balance_temperature_departures()
        dTs[np.abs(dTs) > 15] = 15
        for stage, dT in zip(self.stages, dTs):
            partition = stage.partition
            partition.T += dT
            for i in partition.outs: i.T += dT
       
    def run_mass_balance(self):
        partitions = self.partitions
        Sb, safe = bottoms_stripping_factors_safe(
            np.array([i.B for i in partitions]), 
            np.array([i.K for i in partitions]),
        )
        return top_flow_rates(Sb, *self._iter_args, safe)
       
    def update_mass_balance(self):
        self.set_flow_rates(self.run_mass_balance())
        
    def interpolate_missing_variables(self):
        stages = self.stages
        lle = self._has_lle and 'K' not in self.partition_data
        partitions = [i.partition for i in stages]
        Bs = []
        Ks = []
        Ts = []
        if lle: gamma_y = []
        N_stages = self.N_stages
        index = []
        for i in range(N_stages):
            partition = partitions[i]
            B = partition.B
            T = partition.T
            K = partition.K
            if B is None or K is None: continue
            index.append(i)
            Bs.append(B)
            Ks.append(K)
            Ts.append(T)
            if lle: gamma_y.append(partition.gamma_y)
        N_ok = len(index)
        if len(index) == N_stages:
            Bs = np.array(Bs)
            Ks = np.array(Ks)
            Ts = np.array(Ts)
            if lle: gamma_y = np.array(gamma_y)
        else:
            if N_ok > 1:
                neighbors = get_neighbors(index, size=N_stages)
                try:
                    Bs = fillmissing(neighbors, expand(Bs, index, N_stages))
                except:
                    breakpoint()
                Ts = fillmissing(neighbors, expand(Ts, index, N_stages))
                N_chemicals = self._N_chemicals
                all_Ks = np.zeros([N_stages, N_chemicals])
                if lle: all_gamma_y = all_Ks.copy()
                for i in range(N_chemicals):
                    all_Ks[:, i] = fillmissing(
                        neighbors, 
                        expand([stage[i] for stage in Ks], index, N_stages)
                    )
                    if not lle: continue
                    all_gamma_y[:, i] = fillmissing(
                        neighbors, 
                        expand([stage[i] for stage in gamma_y], index, N_stages)
                    )
                if lle: gamma_y = all_gamma_y
                Ks = all_Ks
            elif N_ok == 1:
                Bs = np.array(N_stages * Bs)
                Ks = np.array(N_stages * Ks)
                Ts = np.array(N_stages * Ts)
                if lle: gamma_y = np.array(N_stages * gamma_y)
            elif N_ok == 0:
                raise RuntimeError('no phase equilibrium')
            for i, stage in enumerate(stages): 
                partition = stage.partition
                T = Ts[i]
                partition.T = T 
                for j in partition.outs: j.T = T
                if partition.B_specification is None: partition.B = Bs[i]
                partition.K = Ks[i]
                if lle: partition.gamma_y = gamma_y[i]
    
    def set_KTBs(self, KTBs):
        if (KTBs <= 0.).any() or not np.isfinite(KTBs).all():
            raise RuntimeError('infeasible equilibrium condition')
        lle = self._has_lle
        N_stages = self.N_stages 
        N_chemicals = self._N_chemicals
        N_flows = N_stages * N_chemicals
        K = KTBs[:N_flows]
        if lle: Ts = KTBs[N_flows:-N_stages]
        Bs = KTBs[-N_stages:]
        K = K.reshape([N_stages, N_chemicals])
        partitions = self.partitions
        N_chemicals = self._N_chemicals
        for i, partition in enumerate(partitions):
            if partition.B_specification is None: partition.B = Bs[i]
            if lle: partition.T = Ts[i]
            partition.K = K[i]
    
    def get_KTBs(self):
        lle = self._has_lle
        N_stages = self.N_stages
        N_chemicals = self._N_chemicals
        N_flows = N_stages * N_chemicals
        KTBs = np.zeros(N_flows + (1 + lle) * N_stages)
        if lle: Ts = KTBs[N_flows:-N_stages]
        Bs = KTBs[-N_stages:]
        last_index = 0
        new_index = N_chemicals
        for i, partition in enumerate(self.partitions):
            KTBs[last_index: new_index] = partition.K
            if lle: Ts[i] = partition.T
            if partition.B_specification is None:
                Bs[i] = partition.B
            else:
                Bs[i] = 0 # It doesnt matter, but it cannot be infinite
            last_index = new_index
            new_index += N_chemicals
        return KTBs
    
    def _overall_error(self, top_flow_rates):
        self._iter(
            top_flow_rates.reshape([self.N_stages, self._N_chemicals])
        ).flatten()
        H_out = np.array([i.H_out for i in self.stages])
        H_in = np.array([i.H_in for i in self.stages])
        diff = H_out - H_in
        diff_mask = np.abs(diff) > 1e-12
        diff = diff[diff_mask]
        denominator = H_out[diff_mask]
        H_in = H_in[diff_mask]
        denominator_mask = np.abs(denominator) < 1e-12
        denominator[denominator_mask] = H_in[denominator_mask]
        errors = diff / denominator
        MSE = (errors * errors).sum()
        return MSE
    
    def _iter(self, variables, KTBs=False):
        self.iter += 1
        if KTBs:
            self.set_KTBs(variables)
            top_flow_rates = self.run_mass_balance()
        else:
            top_flow_rates = variables
        stages = self.stages
        P = self.P
        if self._has_vle:
            self.set_flow_rates(top_flow_rates)
            for i in stages:
                mixer = i.mixer
                partition = i.partition
                mixer.outs[0].mix_from(
                    mixer.ins, energy_balance=False,
                )
                partition._run_decoupled_KTvle(P=P)
                T = partition.T
                for i in (partition.outs + i.outs): i.T = T
            self.update_energy_balance_phase_ratios()
        elif self._has_lle: # LLE
            if 'K' in self.partition_data:
                self.set_flow_rates(top_flow_rates)
            else:
                def psuedo_equilibrium(top_flow_rates):
                        self.set_flow_rates(top_flow_rates)
                        for n, i in enumerate(stages): 
                            mixer = i.mixer
                            partition = i.partition
                            mixer.outs[0].mix_from(
                                mixer.ins, energy_balance=False,
                            )
                            partition._run_decoupled_Kgamma(P=P)
                        return self.run_mass_balance()
                options = tmo.equilibrium.LLE.pseudo_equilibrium_inner_loop_options.copy()
                options['xtol'] *= self.feed_flows.max()
                self.set_flow_rates(
                    flx.fixed_point(
                        psuedo_equilibrium, top_flow_rates, **options,
                    )
                )
            for i in stages: 
                mixer = i.mixer
                partition = i.partition
                mixer.outs[0].mix_from(
                    mixer.ins, energy_balance=False,
                )
                partition._run_decoupled_B()
            self.update_energy_balance_temperatures()
        if self.inside_out and self._has_vle:
            self.update_mass_balance()
            N_stages = self.N_stages
            N_chemicals = self._N_chemicals
            T = np.zeros(N_stages)
            hv = T.copy()
            hl = T.copy()
            specification_index = []
            for i, j in enumerate(self.partitions):
                top, bottom = j.outs
                T[i] = j.T
                if bottom.isempty():
                    top.phase = 'l'
                    hl[i] = top.h
                    top.phase = 'g'
                else:
                    hl[i] = bottom.h
                if top.isempty():
                    bottom.phase = 'g'
                    hv[i] = bottom.h
                    bottom.phase = 'l'
                else:
                    hv[i] = top.h
                if j.B_specification is not None: specification_index.append(i)
            variables = solve_inside_loop(
                self.get_KTBs(), T, hv, hl, self.feed_flows,
                self._asplit_1, self._bsplit_1, 
                self._asplit_left, self._bsplit_left,
                N_stages, np.array(specification_index, int),
                N_chemicals,
                self.feed_enthalpies
            )
            if KTBs:
                return variables
            else:
                self.set_KTBs(variables)
                return self.run_mass_balance()
        elif KTBs:
            return self.get_KTBs()
        else:
            return self.run_mass_balance()

    def get_top_flow_rates_flat(self):
        N_chemicals = self._N_chemicals
        top_flow_rates = np.zeros(self.N_stages * N_chemicals)
        last_index = 0
        new_index = N_chemicals
        partition_index = self._update_index
        for i, partition in enumerate(self.partitions):
            top_flow_rates[last_index: new_index] = partition.outs[0].mol[partition_index]
            last_index = new_index
            new_index = last_index + N_chemicals
        return top_flow_rates
    
    def get_top_flow_rates(self):
        top_flow_rates = np.zeros([self.N_stages, self._N_chemicals])
        partition_index = self._update_index
        for i, partition in enumerate(self.partitions):
            top_flow_rates[i] = partition.outs[0].mol[partition_index]
        return top_flow_rates

    def _conditional_iter(self, top_flow_rates):
        mol = top_flow_rates.flatten()
        top_flow_rates_new = self._iter(top_flow_rates)
        mol_new = top_flow_rates_new.flatten()
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
        return top_flow_rates_new, not_converged

    def _root_iter(self, KTBs):
        KTBs_new = self._iter(
            KTBs, True,
        )
        return KTBs_new - KTBs

    def _sequential_iter(self, top_flow_rates):
        self.fallback_iter += 1
        self.set_flow_rates(top_flow_rates)
        for i in self.stages: i._run()
        for i in reversed(self.stages): i._run()
        mol = top_flow_rates.flatten()
        top_flow_rates = self.get_top_flow_rates()
        mol_new = top_flow_rates.flatten()
        mol_errors = abs(mol - mol_new)
        if mol_errors.any():
            mol_error = mol_errors.max()
            if mol_error > 1e-12:
                nonzero_index, = (mol_errors > 1e-12).nonzero()
                mol_errors = mol_errors[nonzero_index]
                max_errors = np.maximum.reduce([abs(mol[nonzero_index]), abs(mol_new[nonzero_index])])
                rmol_error = (mol_errors / max_errors).max()
                not_converged = (
                    self.fallback_iter < self.fallback_maxiter and (mol_error > self.molar_tolerance
                     or rmol_error > self.relative_molar_tolerance)
                )
            else:
                not_converged = False
        else:
            not_converged = False
        return top_flow_rates, not_converged


# %% General functional algorithms based on MESH equations to solve multi-stage 

@njit(cache=True)
def solve_TDMA(a, b, c, d): # Tridiagonal matrix solver
    """
    Solve a tridiagonal matrix using Thomas' algorithm.
    
    http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    
    Notes
    -----
    `a` array starts from a1 (not a0).
    
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
def solve_TDMA_2D_careful(a, b, c, d, ab_fallback):
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        ai = a[i]
        bi = b[i]
        m = bi.copy()
        inf_mask = bi == inf
        zero_mask = bi == 0
        ok_mask = ~inf_mask & ~zero_mask
        ok_index, = np.nonzero(ok_mask)
        inf_index, = np.nonzero(inf_mask)
        zero_index, = np.nonzero(zero_mask)
        special_index, = np.nonzero(inf_mask & (ai == -inf))
        special = ab_fallback[i]
        for j in inf_index: m[j] = 0
        for j in special_index: m[j] = special
        for j in zero_index: m[j] = inf
        for j in ok_index: m[j] = ai[j] / bi[j]
        b[inext] = b[inext] - m * c[i] 
        d[inext] = d[inext] - m * d[i]
        
    bn = d[n] / b[n]
    bn[bn < 0] = 0
    b[n] = bn
    for i in range(n-1, -1, -1):
        bi = (d[i] - c[i] * b[i+1]) / b[i]
        bi[bi < 0] = 0
        b[i] = bi
    return b

@njit(cache=True)
def solve_LBDMA(a, b, d): # Left bidiagonal matrix solver
    """
    Solve a left bidiagonal matrix using a reformulation of Thomas' algorithm.
    """
    n = d.shape[0] - 1 # number of equations minus 1
    for i in range(n):
        inext = i + 1
        m = a[i] / b[i]
        d[inext] = d[inext] - m * d[i]
    
    b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        b[i] = d[i] / b[i]
    return b

@njit(cache=True)
def solve_RBDMA_1D_careful(b, c, d):
    n = d.shape[0] - 1 # number of equations minus 1
    bn = b[n]
    dn = d[n]
    if bn == 0:
        if dn == 0:
            b[n] = 0
        else:
            b[n] = inf
    else:
        b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        bi = b[i]
        num = d[i] - c[i] * b[i+1]
        if bi == 0:
            if num == 0:
                b[i] = 0
            else:
                b[i] = inf
        else:
            b[i] = num / bi
    return b

@njit(cache=True)
def solve_RBDMA(b, c, d): # Right bidiagonal matrix solver
    """
    Solve a right bidiagonal matrix using a reformulation of Thomas' algorithm.
    """
    n = d.shape[0] - 1 # number of equations minus 1
    b[n] = d[n] / b[n]

    for i in range(n-1, -1, -1):
        b[i] = (d[i] - c[i] * b[i+1]) / b[i]
    return b

@njit(cache=True)
def hot_start_top_flow_rates(
        bottom_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, asplit_1, bsplit_1, N_stages,
    ):
    """
    Solve a-phase flow rates for a single component across 
    equilibrium stages with side draws. 

    Parameters
    ----------
    bottom_flows : Iterable[1d array]
        Bottom flow rates by stages.
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a (extract or vapor) and b (raffinate or liquid) leaving the stage 
        respectively.
    stage_index : 1d array
        Stage index for phase ratios.
    top_feed_flows : Iterable[1d array]
        Top flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    bottom_feed_flows : Iterable [1d array]
        Bottom flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = top_feed_flows.copy()
    b = d.copy()
    c = d.copy()
    for i in range(N_stages): 
        c[i] = asplit_1[i]
        b[i] = 1
    for n in range(stage_index.size):
        i = stage_index[n]
        B = phase_ratios[n]
        if B <= 1e-32:
            b[i] = inf
        else:
            b[i] += 1 / B 
        if i == 0:
            d[i] += bottom_feed_flows[i]
        else:
            d[i] += bottom_feed_flows[i] - bottom_flows[i - 1] * bsplit_1[i - 1]
    return solve_RBDMA(b, c, d)

@njit(cache=True)
def hot_start_bottom_flow_rates(
        top_flows, phase_ratios, stage_index, top_feed_flows,
        bottom_feed_flows, asplit_1, bsplit_1, N_stages
    ):
    """
    Solve a-phase flow rates for a single component across 
    equilibrium stages with side draws. 

    Parameters
    ----------
    bottom_flows : Iterable[1d array]
        Bottom flow rates by stages.
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a (extract or vapor) and b (raffinate or liquid) leaving the stage 
        respectively.
    stage_index : 1d array
        Stage index for phase ratios.
    top_feed_flows : Iterable[1d array]
        Top flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    bottom_feed_flows : Iterable [1d array]
        Bottom flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a: 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    d = bottom_feed_flows.copy()
    b = d.copy()
    a = d.copy()
    for i in range(N_stages): 
        a[i] = bsplit_1[i]
        b[i] = 1
    last_stage = N_stages - 1
    for n in range(stage_index.size):
        i = stage_index[n]
        b[i] += phase_ratios[n]
        if i == last_stage:
            d[i] += top_feed_flows[i]
        else:
            d[i] += top_feed_flows[i] - top_flows[i + 1] * asplit_1[i + 1]
    return solve_LBDMA(a, b, d)

@njit(cache=True)
def bottoms_stripping_factors_safe(phase_ratios, partition_coefficients):
    """
    Return the bottoms stripping factors (i.e., the ratio of components in 
    the bottoms over the top) and a flag dictating whether it is safe for division
    and multiplication (i.e., whether 0 or inf are present).
    
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

    """
    zero_mask = phase_ratios <= 0.
    inf_mask = phase_ratios >= 1e32
    ok_mask = ~zero_mask & ~inf_mask
    phase_ratios = np.expand_dims(phase_ratios, -1)
    safe = ok_mask.all()
    if safe:
        # Bottoms stripping factor are, by definition, the ratio of components in the bottoms over the top.
        bottoms_stripping_factors = 1. / (phase_ratios * partition_coefficients)
    else:
        zero_index, = np.nonzero(zero_mask)
        inf_index, = np.nonzero(inf_mask)
        ok_index, = np.nonzero(ok_mask)
        bottoms_stripping_factors = np.zeros(partition_coefficients.shape)
        for i in ok_index:
            bottoms_stripping_factors[i] = 1. / (phase_ratios[i] * partition_coefficients[i])
        for i in zero_index:
            bottoms_stripping_factors[i] = inf
        for i in inf_index:
            bottoms_stripping_factors[i] = 0.
    return bottoms_stripping_factors, safe

@njit(cache=True)
def top_flow_rates(
        bottoms_stripping_factors, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
        safe,
    ):
    """
    Solve a-phase flow rates for a single component across equilibrium stages with side draws. 

    Parameters
    ----------
    bottoms_stripping_factors : Iterable[1d array]
        The ratio of component flow rates in phase b (raffinate or liquid) over
        the component flow rates in phase a (extract or vapor). 
    feed_flows : Iterable[1d array]
        Flow rates of all components fed across stages. Shape should be 
        (N_stages, N_chemicals).
    asplit_1 : 1d array
        Side draw split from phase a minus 1 by stage.
    bsplit_1 : 1d array
        Side draw split from phase b minus 1 by stage.

    Returns
    -------
    flow_rates_a : 2d array
        Flow rates of phase a with stages by row and components by column.

    """
    b = 1. + bottoms_stripping_factors
    c = asplit_1[1:]
    d = feed_flows.copy()
    a = np.expand_dims(bsplit_1, -1) * bottoms_stripping_factors
    if safe:    
        top_flows = solve_TDMA(a, b, c, d) 
    else:
        top_flows = solve_TDMA_2D_careful(a, b, c, d, bsplit_1)
    return top_flows

@njit(cache=True)
def bottom_flow_rates(
        top_flows, feed_flows, asplit_left, bsplit_left,
        N_stages, N_chemicals
    ):
    bottom_flows = 0 * top_flows
    bottom_flows[0] = feed_flows[0] + top_flows[1] * asplit_left[1] - top_flows[0]
    for i in range(1, N_stages-1):
        bottom_flows[i] = (
            feed_flows[i] + bsplit_left[i-1] * bottom_flows[i-1] + 
            top_flows[i+1] * asplit_left[i+1] - top_flows[i]
        )
    bottom_flows[-1] = feed_flows[-1] + bsplit_left[-2] * bottom_flows[-2] - top_flows[-1]
    return bottom_flows

@njit(cache=True)
def mass_balance(
        top_flows, feed_flows, asplit_left, bsplit_left,
        correct_stages, N_stages, N_chemicals
    ):
    bottom_flows = 0 * top_flows
    infeasible = True
    while infeasible:
        row = feed_flows[0] + top_flows[1] * asplit_left[1] - top_flows[0]
        index = row < 0
        infeasible = index.any()
        if infeasible:
            infeasible_flows = row[index]
            row[index] = 0
            if correct_stages[0]:
                top_flows[0][index] += infeasible_flows 
                correct_stages[0] = False
                break
            else:
                infeasible = False
        bottom_flows[0] = row
        for i in range(1, N_stages-1):
            row = (
                feed_flows[i] + bsplit_left[i-1] * bottom_flows[i-1] + 
                top_flows[i+1] * asplit_left[i+1] - top_flows[i]
            )
            index = row < 0
            infeasible = index.any()
            if infeasible:
                infeasible_flows = row[index]
                row[index] = 0
                if correct_stages[i]:
                    top_flows[i][index] += infeasible_flows 
                    correct_stages[i] = False
                    break
                else:
                    infeasible = False
            bottom_flows[i] = row
        row = feed_flows[-1] + bsplit_left[-2] * bottom_flows[-2] - top_flows[-1]
        index = row < 0
        infeasible = index.any()
        if infeasible:
            infeasible_flows = row[index]
            row[index] = 0
            if correct_stages[-1]:
                top_flows[-1][index] += infeasible_flows 
                correct_stages[-1] = False
                break
            else:
                infeasible = False
        bottom_flows[-1] = row
    return bottom_flows

# @njit(cache=True)
# def balance_flow_rates(top_flows, bottom_flows, asplit_left, bsplit_left):
#     pass
    
@njit(cache=True)
def phase_ratio_departures(
        L, V, hl, hv, asplit_1, asplit_left, bsplit_left, 
        N_stages, specification_index, H_feeds
    ):
    # hV1*L1*dB1 - hv2*L2*dB2 = Q1 + H_in - H_out
    b = hv * L
    c = b[1:] * asplit_1[1:]
    Hl_out = hl * L
    Hv_out = hv * V
    d = H_feeds - Hl_out - Hv_out
    Hl_in = (Hl_out * bsplit_left)[:-1]
    Hv_in = (Hv_out * asplit_left)[1:]
    d[1:] += Hl_in
    d[:-1] += Hv_in
    for i, j in enumerate(specification_index):
        b[j] = 0
        d[j] = 0
        jlast = j - 1
        if jlast > 0: c[jlast] = 0
        try: c[j] = 0
        except: pass
    return solve_RBDMA_1D_careful(b, c, d)

@njit(cache=True)
def temperature_departures(Cv, Cl, Hv, Hl, asplit_left, bsplit_left,
                           N_stages, H_feeds):
    # ENERGY BALANCE
    # C1dT1 - Cv2*dT2 - Cl0*dT0 = Q1 - H_out + H_in
    b = (Cv + Cl)
    a = -(Cl * bsplit_left)
    c = -(Cv * asplit_left)[1:]
    d = H_feeds - Hl - Hv
    d[1:] += (Hl * bsplit_left)[:-1]
    d[:-1] += (Hv * asplit_left)[1:]
    return solve_TDMA(a, b, c, d)

def get_neighbors(index=None, all_index=None, missing=None, size=None):
    if size is not None: all_index = set(range(size))
    if sum([i is None for i in (index, all_index, missing)]) > 1:
        raise ValueError('at least two arguments must be given')
    if missing is None: 
        missing = set(all_index).difference(index)
    else:
        missing = set(missing)
    if index is None:
        index = all_index.difference(missing)
    else:
        index = set(index)
    if all_index is None:
        all_index = index | missing
    else:
        all_index = set(all_index)
    size = len(all_index)
    index_set = set(index)
    neighbors = []
    for i in missing:
        lb = i
        while lb > -1:
            lb -= 1
            if lb in index_set: break
        ub = i
        while ub < size:
            ub += 1
            if ub in index_set: break
        if ub == size:
            neighbors.append(
                (i, (lb,))
            )
        elif lb == -1:
            neighbors.append(
                (i, (ub,))
            )
        else:
            neighbors.append(
                (i, (lb, ub))
            )
    return neighbors

def expand(values, index, size):
    new_values = np.zeros(size)
    new_values[index] = values
    return new_values

def fillmissing(all_neighbors, values):
    for i, neighbors in all_neighbors:
        if len(neighbors) == 2:
            lb, ub = neighbors
            lb_distance = i - lb
            ub_distance = ub - i
            sum_distance = lb_distance + ub_distance
            wlb = ub_distance / sum_distance
            wub = lb_distance / sum_distance
            x = wlb * values[lb] + wub * values[ub]
            values[i] = x
        else:
            values[i] = values[neighbors[0]]
    return values

# %% Methods for root finding

options = dict(ftol=1e-3, maxiter=100)
for name in ('anderson', 'diagbroyden', 'excitingmixing', 'linearmixing', 
             'broyden1', 'broyden2', 'krylov', 'hybr'):
    MultiStageEquilibrium.root_options[name] = (root, False, {'method': name, 'options': options})

# %% Russel's inside-out algorithm

@njit(cache=True)
def omega_approx(y, K):
    y_over_K = (y / K)
    return y_over_K / y_over_K.sum()

@njit(cache=True)
def Kb_init(y, K):
    omega = omega_approx(y, K)
    return np.exp((omega * np.log(K)).sum(axis=1))

@njit(cache=True)
def Kb_iter(alpha, x):
    return 1 / (alpha * x).sum(axis=1)

@njit(cache=True)
def alpha_approx(K, Kb):
    return K / Kb

@njit(cache=True)
def fit(x, y):
    xmean = x.mean()
    ymean = y.mean()
    xxmean = x - xmean
    m = (xxmean * (y - ymean)).sum() / (xxmean * xxmean).sum()
    b = ymean - m * xmean
    return m, b

@njit(cache=True)
def fit_partition_model(T, Kb):
    x = 1 / T
    y = np.log(Kb)
    xdiff = np.diff(x)
    ydiff = np.diff(y)
    M = y.copy()
    M[:-1] = ydiff / xdiff
    M[-1] = M[-2]
    M[-2] = (M[-3] + M[-1]) / 2
    B = y - M * x
    return M, B

@njit(cache=True)
def h_approx(T, m, b):
    return m * T + b

@njit(cache=True)
def T_approx(Kb, m, b):
    return m / (np.log(Kb) - b)

def solve_inside_loop(
        KB, T, hv, hl, feed_flows,
        asplit_1, bsplit_1, asplit_left, bsplit_left,
        N_stages, specification_index, N_chemicals, H_feeds
    ):
    correct_stages = np.ones(N_stages, bool)
    args = inside_loop_args(
            KB, T, hv, hl, feed_flows,
            asplit_1, bsplit_1, asplit_left, bsplit_left,
            N_stages, specification_index, N_chemicals, H_feeds,
            correct_stages
        )
    # result = root(inside_loop, KB.flatten(), 
    #               options=dict(ftol=1e-6), args=args)
    # print(result.x)
    KB_new = flx.fixed_point(inside_loop, KB.flatten(), xtol=1e-6, args=args, checkiter=False)
    # print(KB_new)
    return KB_new
   
@njit(cache=True)
def inside_loop_args(
        KB, T, hv, hl, feed_flows,
        asplit_1, bsplit_1, asplit_left, bsplit_left,
        N_stages, specification_index, N_chemicals, H_feeds,
        correct_stages,
    ):
    N_flows = N_stages * N_chemicals
    K = KB[:N_flows]
    B = KB[N_flows:]
    K = K.reshape((N_stages, N_chemicals))
    Sb, safe = bottoms_stripping_factors_safe(B, K)
    top_flows = top_flow_rates(
        Sb, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
        safe,
    )
    dummy = top_flows.sum(axis=1)
    dummy[dummy == 0] = 1
    y = top_flows / np.expand_dims(dummy, -1)
    Kb = Kb_init(y, K)
    Kb_coef = fit_partition_model(T, Kb)
    hv_coef = fit(T, hv)
    hl_coef = fit(T, hl)
    alpha = alpha_approx(K, np.expand_dims(Kb, -1))
    return (alpha, Kb_coef, hv_coef, hl_coef, 
            feed_flows, asplit_1, bsplit_1,
            asplit_left, bsplit_left, N_stages, 
            specification_index, N_chemicals, H_feeds,
            correct_stages)
    
@njit(cache=True)
def inside_loop(KB, alpha, Kb_coef, hv_coef, hl_coef, 
                feed_flows, asplit_1, bsplit_1,
                asplit_left, bsplit_left, N_stages, 
                specification_index, N_chemicals, H_feeds,
                correct_stages):
    N_flows = N_stages * N_chemicals
    K = KB[:N_flows]
    B = KB[N_flows:]
    K = K.reshape((N_stages, N_chemicals))
    Sb, safe = bottoms_stripping_factors_safe(B, K)
    top_flows = top_flow_rates(
        Sb, 
        feed_flows,
        asplit_1,
        bsplit_1,
        N_stages,
        safe,
    )
    bottom_flows = mass_balance(
        top_flows, feed_flows, asplit_left, bsplit_left, 
        correct_stages.copy(), N_stages, N_chemicals
    )
    top_flows_net = top_flows.sum(axis=1)
    bottom_flows_net = bottom_flows.sum(axis=1)
    dummy = bottom_flows_net.copy()
    dummy[dummy == 0] = 1e-12
    x = bottom_flows / dummy[:, np.newaxis]
    dummy = top_flows_net.copy()
    dummy[dummy == 0] = 1e-12
    y = top_flows / dummy[:, np.newaxis]
    mask = (x == 0)
    for i in mask: x[i] =  y[i] / K[i]
    Kb = Kb_iter(alpha, x)
    KB_new = KB.copy()
    last_index = 0
    new_index = N_chemicals
    for row in (alpha * Kb[:, np.newaxis]):
        KB_new[last_index: new_index] = row
        new_index = last_index + N_chemicals
    T = T_approx(Kb, *Kb_coef)
    T.sort()
    hv = h_approx(T, *hv_coef)
    hl = h_approx(T, *hl_coef)
    KB_new[N_flows:] = B + phase_ratio_departures(
        bottom_flows_net, top_flows_net, hl, hv, asplit_1, 
        asplit_left, bsplit_left, N_stages,
        specification_index, H_feeds
    )
    return KB_new
    
    
    
    