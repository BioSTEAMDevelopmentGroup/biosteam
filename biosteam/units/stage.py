# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains abstract classes for modeling stage-wise separations/reactions in unit operations.

"""
import thermosteam as tmo
from thermosteam.utils import jitdata
from thermosteam.base.sparse import SparseVector, sum_sparse_vectors
from thermosteam import separations as sep
from numba import njit
import biosteam as bst
import flexsolve as flx
import numpy as np
import pandas as pd
from scipy.optimize import minimize, differential_evolution
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator, PchipInterpolator, Akima1DInterpolator
from math import inf
from typing import Callable, Optional
from copy import copy
from scipy.optimize import fsolve
from scipy.optimize._numdiff import approx_derivative
from scipy.differentiate import jacobian
from .. import Unit
from .design_tools import MESH
from numba import float64, types, njit
from thermosteam import (
    equilibrium, VariableNode,
)

__all__ = (
    'SinglePhaseStage',
    'ReactivePhaseStage',
    'StageEquilibrium',
    'MultiStageEquilibrium',
    'PhasePartition',
)

# %% Equation-oriented tools

@jitdata
class JacobianData:
    dEdx: float64[:, :]
    dHdFtop: float64[:]
    dHdFbot: float64[:]
    dHdTtop: float
    dHdTbot: float
    split_top: float
    split_bot: float
    variable: str
    value: float


@jitdata
class InterstageData:
    H: float
    mol: float64[:]
    split: float
    
    
@jitdata
class FeedData:
    H: float
    mol: float64[:]
    

@jitdata
class StageData:
    T: float
    top: InterstageData
    bottom: InterstageData
    feed: FeedData

# We can order the variables following standard convention.
# Variable order for reference stage j given n chemicals:
# [A] j-1 (previous stage), [B] j (reference stage), [C] j+1 (bottom stage) 
# For each A, B, and C:
# s_top_i, ... s_top_n, T, s_bottom_i, ... s_bottom_n.

@jitdata
class EquationIndex:
    H: int
    M: types.slice2_type
    E: types.slice2_type
    
    
@jitdata
class VariableIndex:
    Ftop: types.slice2_type
    Fbot: types.slice2_type
    T: int


@jitdata
class JacobianConstructor:
    equation_index: EquationIndex
    variable_index: VariableIndex
    
    def __init__(self, N_chemicals: int):
        N_plus_1 = N_chemicals + 1
        end = N_plus_1 + N_chemicals
        self.equation_index = EquationIndex(
            0,
            slice(1, N_plus_1),
            slice(N_plus_1, end),
        )
        self.variable_index = VariableIndex(
            slice(N_chemicals),
            slice(N_plus_1, end),
            N_chemicals,
        )
    
    def fill_A(self, 
            A: float64[:, :], 
            upper: JacobianData, 
            center: JacobianData
        ):
        variable = center.variable
        eq = self.equation_index
        var = self.variable_index
        if upper.split_bot != 0:
            split = upper.split_bot - 1
            if variable == 'Q':
                # Otherwise, zeros, energy balance is decoupled from top stage
                A[eq.H, var.Fbot] = split * upper.dHdFbot
                A[eq.H, var.T] = split * upper.dHdTbot
            np.fill_diagonal(A[eq.M, var.Fbot], split)
        else:
            if variable == 'Q':
                # Otherwise, zeros, energy balance is decoupled from top stage
                A[eq.H, var.Fbot] = -upper.dHdFbot 
                A[eq.H, var.T] = -upper.dHdTbot
            np.fill_diagonal(A[eq.M, var.Fbot], -1)
        
    def fill_B(self, 
            B: float64[:, :], 
            center: JacobianData
        ):
        variable = center.variable
        eq = self.equation_index
        var = self.variable_index
        if variable == 'Q':
            B[eq.H, var.Ftop] = center.dHdFtop
            B[eq.H, var.Fbot] = center.dHdFbot
            B[eq.H, var.T] = center.dHdTtop + center.dHdTbot
        elif variable == 'T':
            B[eq.H, var.T] = 1
        elif variable == 'B':
            B[eq.H, var.Ftop] = 1
            B[eq.H, var.Fbot] = -center.value
        else:
            raise ValueError("invalid specification variable '" + variable + "'")
        np.fill_diagonal(B[eq.M, var.Ftop], 1)
        np.fill_diagonal(B[eq.M, var.Fbot], 1)
        B[eq.E] = center.dEdx
    
    def fill_C(self, 
            C: float64[:, :], 
            center: JacobianData, 
            lower: JacobianData
        ):
        variable = center.variable
        eq = self.equation_index
        var = self.variable_index
        if lower.split_top:
            split = lower.split_top - 1
            if variable == 'Q':
                # Otherwise, zeros, energy balance is decoupled from bottom stage
                C[eq.H, var.Ftop] = split * lower.dHdFtop
                C[eq.H, var.T] = split * lower.dHdTtop
            np.fill_diagonal(C[eq.M, var.Ftop], split)
        else:
            if variable == 'Q':
                # Otherwise, zeros, energy balance is decoupled from bottom stage
                C[eq.H, var.Ftop] = -lower.dHdFtop
                C[eq.H, var.T] = -lower.dHdTtop
            np.fill_diagonal(C[eq.M, var.Ftop], -1)

# @njit(cache=True)
def jacobian_blocks(jacobian_data, N_stages, N_chemicals, N_variables):
    JC = JacobianConstructor(N_chemicals)
    A_blocks = np.zeros((N_stages-1, N_variables, N_variables))
    B_blocks = np.zeros((N_stages, N_variables, N_variables))
    C_blocks = np.zeros((N_stages-1, N_variables, N_variables))
    center = jacobian_data[0]
    lower = jacobian_data[1]
    JC.fill_B(B_blocks[0], center)
    JC.fill_C(C_blocks[0], center, lower)
    for i in range(1, N_stages-1):
        upper = center
        center = lower
        lower = jacobian_data[i]
        JC.fill_A(A_blocks[i], upper, center)
        JC.fill_B(B_blocks[i], center)
        JC.fill_C(C_blocks[i], center, lower)
    upper = center
    center = lower
    JC.fill_A(A_blocks[-1], upper, center)
    JC.fill_B(B_blocks[-1], center)
    return A_blocks, B_blocks, C_blocks


# %% Single phase, only partially supported in BioSTEAM (no active tests)

class SinglePhaseStage(Unit):
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    
    @property
    def _energy_variable(self):
        if self.T is None: return 'T'
        else: return None
    
    def _init(self, T=None, P=None, Q=None, phase=None):
        self.T = T
        self.Q = Q
        self.P = P
        self.phase = phase
        
    def _run(self):
        outlet = self.outs[0]
        outlet.mix_from(self.ins, energy_balance=False)
        if self.P is not None: outlet.P = self.P
        if self.phase is None: 
            outlet.phase = self.ins[0].phase
        else:
            outlet.phase = self.phase
        if self.T is None:
            if self.Q is None:
                raise RuntimeError('must specify either Q or T')
            else:
                outlet.H = sum([i.H for i in self.ins], self.Q)
        elif self.Q is not None:
            raise RuntimeError('must specify either Q or T; not both')
        else:
            outlet.T = self.T

    def _get_energy_departure_coefficient(self, stream):
        if self.T is None: return (self, -stream.C)
    
    def _create_energy_departure_equations(self):
        if self.T is not None: return []
        # Ll: C1dT1 - Ce2*dT2 - Cr0*dT0 - hv2*L2*dB2 = Q1 - H_out + H_in
        # gl: hV1*L1*dB1 - hv2*L2*dB2 - Ce2*dT2 - Cr0*dT0 = Q1 + H_in - H_out
        outlet = self.outs[0]
        coeff = {self: outlet.C}
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        return [(coeff, outlet.H - sum([i.H for i in self.ins]))]
        
    def _create_material_balance_equations(self, composition_sensitive):
        outlet = self.outs[0]
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {outlet: ones}
        for i in process_inlets: 
            if i in eq_overall:
                del eq_overall[i]
            else:
                eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        )
        return equations

    def _update_energy_variable(self, departure):
        for i in self.outs: i.T += departure
        
    def _update_nonlinearities(self): pass
    
    @property
    def equation_node_names(self): 
        if self._energy_variable is None:
            return (
                'overall_material_balance_node', 
            )
        else:
            return (
                'overall_material_balance_node', 
                'energy_balance_node',
            )
    
    def initialize_overall_material_balance_node(self):
        self.overall_material_balance_node.set_equations(
            outputs=[i.F_node for i in self.outs],
            inputs=[j for i in self.ins if (j:=i.F_node)],
        )
        
    def initialize_energy_balance_node(self):
        self.energy_balance_node.set_equations(
            inputs=(
                self.T_node, 
                *[i.T_node for i in (*self.ins, *self.outs)],
                *[i.F_node for i in (*self.ins, *self.outs)],
                *[j for i in self.ins if (j:=i.E_node)],
            ),
            outputs=[
                j for i in self.outs if (j:=i.E_node)
            ],
        )
    
    @property
    def T_node(self):
        if hasattr(self, '_T_node'): return self._T_node
        self._T_node = var = VariableNode(f"{self.node_tag}.T", lambda: self.T)
        return var 
        
    @property
    def E_node(self):
        if self._energy_variable is None:
            return None
        else:
            return self.T_node

    def _collect_edge_errors(self):
        equation_name = self.overall_material_balance_node.name
        outs = self.outs
        IDs = self.chemicals.IDs
        results = []
        error = sum([i.mol for i in outs]) - sum([i.mol for i in self.ins])
        for i, outlet in enumerate(outs):
            for j, ID in enumerate(IDs):
                index = (equation_name, outlet.F_node.name, ID)
                results.append((index, error[j]))
        return results # list[tuple[tuple[equation_name, variable_name, chemical_name | '-'], value]]

    def _collect_equation_errors(self):
        equation_name = self.overall_material_balance_node.name
        outs = self.outs
        results = []
        error = np.abs(sum([i.mol for i in outs]) - sum([i.mol for i in self.ins])).sum()
        index = equation_name
        results.append((index, error))
        
        if self._energy_variable is not None:
            equation_name = self.energy_balance_node.name
            error = sum([i.H for i in outs]) - sum([i.H for i in self.ins])
            results.append((equation_name, np.abs(error)))
        
        return results # list[tuple[equation_name, value]]
    

class ReactivePhaseStage(bst.Unit): # Does not include VLE
    _N_outs = _N_ins = 1
    _ins_size_is_fixed = False
    
    @property
    def equation_node_names(self): 
        if self._energy_variable is None:
            return (
                'overall_material_balance_node',
                'reaction_phenomenode',
            )
        else:
            return (
                'overall_material_balance_node', 
                'reaction_phenomenode',
                'energy_balance_node',
            )
    
    @property 
    def _energy_variable(self):
        if self.T is None: return 'T'
    
    def _init(self, reaction, T=None, P=None, Q=0, phase=None):
        self.reaction = reaction
        self.T = T
        self.P = P
        self.Q = Q
        self.phase = phase
        
    def _run(self):
        feed = self.ins[0]
        outlet, = self.outs
        outlet.copy_like(feed)
        if self.P is not None: outlet.P = self.P
        if self.phase is not None: outlet.phase = self.phase
        if self.T is None: 
            self.reaction.adiabatic_reaction(outlet, Q=self.Q)
        else:
            self.reaction(outlet)
            outlet.T = self.T
        self.dmol = outlet.mol - feed.mol
        
    def _create_material_balance_equations(self, composition_sensitive=False):
        product, = self.outs
        n = self.chemicals.size
        ones = np.ones(n)
        minus_ones = -ones
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        # Overall flows
        eq_overall = {}
        predetermined_flow = SparseVector.from_dict(sum_sparse_vectors([i.mol for i in fresh_inlets]), size=n)
        rhs = predetermined_flow + self.dmol
        eq_overall[product] = ones
        for i in process_inlets: eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, rhs)
        )
        return equations
    
    def _update_energy_variable(self, departure):
        self.outs[0].T += departure
    
    def _get_energy_departure_coefficient(self, stream):
        if self.T is not None: return
        return (self, -stream.C)
    
    def _create_energy_departure_equations(self):
        if self.T is not None: return []
        # Ll: C1dT1 - Ce2*dT2 - Cr0*dT0 - hv2*L2*dB2 = Q1 - H_out + H_in
        # gl: hV1*L1*dB1 - hv2*L2*dB2 - Ce2*dT2 - Cr0*dT0 = Q1 + H_in - H_out
        outlet = self.outs[0]
        coeff = {self: outlet.C}
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        return [(coeff, -self.Hnet)]
    
    def _update_nonlinearities(self):
        f = PhasePartition.dmol_relaxation_factor
        old = self.dmol
        new = self.reaction.conversion(self.ins[0])
        self.dmol = f * old + (1 - f) * new
    
    def initialize_reaction_phenomenode(self):
        self.reaction_phenomenode.set_equations(
            inputs=[j for i in self.ins if (j:=i.F_node)],
            outputs=[self.R_node],
        )
    
    def initialize_overall_material_balance_node(self):
        self.overall_material_balance_node.set_equations(
            inputs=[j for i in self.ins if (j:=i.F_node)] + [self.R_node],
            outputs=[i.F_node for i in self.outs],
        )
            
    def initialize_energy_balance_node(self):
        if self.T is None:
            self.energy_balance_node.set_equations(
                inputs=(
                    self.T_node, 
                    *[i.T_node for i in (*self.ins, *self.outs)],
                    *[i.F_node for i in (*self.ins, *self.outs)],
                    *[j for i in self.ins if (j:=i.E_node)]
                ),
                outputs=[j for i in self.outs if (j:=i.E_node)],
            )
    
    @property
    def R_node(self):
        if hasattr(self, '_R_node'): return self._R_node
        self._R_node = var = VariableNode(f"{self.node_tag}.R", lambda: self.dmol)
        return var 
    
    @property
    def T_node(self):
        if hasattr(self, '_T_node'): return self._T_node
        if self.T is None: 
            var = VariableNode(f"{self.node_tag}.T", lambda: self.T)
        else:
            var = None
        self._T_node = var
        return var 
    
    def get_E_node(self, stream):
        if self.T is None:
            return self.E_node
        else:
            return None
    
    @property
    def E_node(self):
        if hasattr(self, '_E_node'): return self._E_node
        if self._energy_variable is None:
            var = None
        else:
            var = self.T_node
        self._E_node = var
        return var


# %% Two phases (fully supported in BioSTEAM)

class StageEquilibrium(Unit):
    _N_ins = 0
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    auxiliary_unit_names = ('partition', 'mixer', 'splitters')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
            phases, partition_data=None, top_split=0, bottom_split=0,
            top_chemical=None, reaction=None, **specifications
        ):
        self._N_outs = 2 + int(top_split) + int(bottom_split)
        self.phases = phases
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
        self.reaction = reaction
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
        self.specify_variables(**specifications)
    
    # %% Streams
    
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
    
    def add_feed(self, stream):
        self.ins.append(stream)
        self.mixer.ins.append(
            self.auxlet(
                stream
            )
        )
    
    # %% Variables and specifications
    
    @property
    def specified_variable(self):
        return self.partition.specified_variable
    @specified_variable.setter
    def specified_variable(self, specified_variable):
        self.partition.specified_variable = specified_variable
    
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
    def T(self):
        return self.partition.T
    @T.setter
    def T(self, T):
        self.partition.T = T
        for i in self.partition.outs: i.T = T
    
    @property
    def P(self):
        return self.partition.P
    @P.setter
    def P(self, P):
        self.partition.P = P
        for i in self.partition.outs: i.P = P
    
    @property
    def K(self):
        return self.partition.K
    @K.setter
    def K(self, K):
        self.partition.K = K
    
    @property
    def reaction(self):
        return self.partition.reaction
    @reaction.setter
    def reaction(self, reaction):
        self.partition.reaction = reaction
        
    def specify_variables(self, **specifications):
        partition = self.partition
        partition.Q = partition.B = partition.T = partition.F = partition.P = None
        specified_variable = None
        P_specified = False
        for name, value in specifications.items():
            if value is None: return
            if name in ('Q', 'Duty'):
                specified_variable = 'Q'
                partition.Q = value
            elif name == 'Reflux':
                partition.B = inf if value == 0 else 1 / value
                specified_variable = 'B'
            elif name in ('B', 'Boilup'):
                partition.B = value
                specified_variable = 'B'
            elif name in ('T', 'Temperature'):
                partition.T = value
                specified_variable = 'T'
            elif name in ('F', 'Flow'):
                partition.F = value
                specified_variable = 'F'
                raise NotImplementedError('F specification not implemented in BioSTEAM yet')
            elif name in ('P', 'Pressure'):
                partition.P = value
                P_specified = True
            else:
                raise RuntimeError(f"specification '{name}' not implemented for stage")
        N_specs = len(specifications)
        if N_specs == 0:
            specified_variable = 'Q'
            self.P = 101325
            self.Q = 0
        elif N_specs == 1:
            if specified_variable is None:
                specified_variable = 'Q'
                self.Q = 0
            else:
                self.P = 101325
        elif N_specs == 2:
            if not P_specified: raise ValueError('can only specify one of Q, T, F, B; 2 specified')
        elif N_specs > 2:
            raise ValueError('cannot specify over 2 variables')
        if specified_variable == 'Q':
            if self.phases == ('g', 'l'):
                self._energy_variable = 'B'
            else:
                self._energy_variable = 'T'
        else:
            self._energy_variable = None
        self.specified_variable = specified_variable
    
    # %% Modular simulation
    
    def _run(self):
        if self.specified_variable == 'T':
            mix = self.mixer.outs[0]
            mix.phase = 'l'
            mix.mol = sum([i.mol for i in self.ins])
            mix.T = self.T
        else:
            self.mixer._run()
        self.partition._run()
        for i in self.splitters: i._run()
        self._update_separation_factors()
    
    def _init_separation_factors(self):
        if self.specified_variable == 'B':
            B = self.B
            if B == 0:
                self.S = 0 * self.K
            elif B == inf:
                self.S = self.K.copy()
                self.S[:] = inf
            else:
                self.S = (B * self.K)
        else:
            B = self.B
            if B is None:
                self.S = np.ones(self.chemicals.size)
            else:
                self.S = (B * self.K)
        
    def _update_separation_factors(self, f=None):
        if not hasattr(self, 'S'): self._init_separation_factors()
        if self.B is None or self.B == inf or self.B == 0: return
        K = self.K
        if K.size != self.S.size: 
            self._init_separation_factors()
            return
        S = K * self.B
        if self.phases == ('L', 'l'):
            self.S = S
        else:
            if f is None: f = self.partition.S_relaxation_factor
            if f == 0:
                self.S = S
            else:
                self.S = np.exp((1 - f) * np.log(S) + f * np.log(self.S)) if f else S
    
    # %% Equation-oriented simulation
    
    def _stage_data(self, x, H_feed, mol_feed, H_magnitude):
        IDs = self.partition.IDs
        chemicals = self.chemicals[IDs]
        P = self.P
        phase_top, phase_bot = self.phases
        N_chemicals = len(IDs)
        mol_top = x[:N_chemicals]
        T = x[N_chemicals]
        mol_bot = x[-N_chemicals:]
        # TODO: Use property package to estimate H
        htop = np.array([i.H(phase_top, T, P) for i in chemicals]) / H_magnitude
        hbot = np.array([i.H(phase_bot, T, P) for i in chemicals]) / H_magnitude
        Htop = (htop * mol_top).sum()
        Hbot = (hbot * mol_bot).sum()
        return StageData(
            T=T,
            top=InterstageData(Htop, mol_top, self.top_split),
            bottom=InterstageData(Hbot, mol_bot, self.bottom_split),
            feed=FeedData(H_feed / H_magnitude, mol_feed),
        )
    
    def _jacobian_data(self, x, H_magnitude):
        IDs = self.partition.IDs
        chemicals = self.chemicals[IDs]
        P = self.P
        phase_top, phase_bot = self.phases
        N_chemicals = len(IDs)
        mol_top = x[:N_chemicals]
        T = x[N_chemicals]
        mol_bot = x[-N_chemicals:]
        Ctop = (np.array([i.Cn(phase_top, T, P) for i in chemicals]) * mol_top).sum() / H_magnitude
        Cbot = (np.array([i.Cn(phase_bot, T, P) for i in chemicals]) * mol_bot).sum() / H_magnitude
        # dEdx = jacobian(self._equilibrium_residuals_vectorized, x).df
        dEdx = approx_derivative(self._equilibrium_residuals_vectorized, x)
        htop = np.array([i.H(phase_top, T, P) for i in chemicals]) / H_magnitude
        hbot = np.array([i.H(phase_bot, T, P) for i in chemicals]) / H_magnitude
        JD = JacobianData(
            dEdx=dEdx,
            dHdFtop=htop,
            dHdFbot=hbot,
            dHdTtop=Ctop,
            dHdTbot=Cbot,
            split_top=self.top_split,
            split_bot=self.bottom_split,
            variable=self.specified_variable,
            value=getattr(self, self.specified_variable),
        )
        return JD
    
    def _simultaneous_correction_data(self, x, H_feed, mol_feed, H_magnitude):
        IDs = self.partition.IDs
        chemicals = self.chemicals[IDs]
        P = self.P
        phase_top, phase_bot = self.phases
        N_chemicals = len(IDs)
        mol_top = x[:N_chemicals]
        T = x[N_chemicals]
        mol_bot = x[-N_chemicals:]
        Ctop = (np.array([i.Cn(phase_top, T, P) for i in chemicals]) * mol_top).sum() / H_magnitude
        Cbot = (np.array([i.Cn(phase_bot, T, P) for i in chemicals]) * mol_bot).sum() / H_magnitude
        # dEdx = jacobian(self._equilibrium_residuals_vectorized, x).df
        dEdx = approx_derivative(self._equilibrium_residuals_vectorized, x)
        htop = np.array([i.H(phase_top, T, P) for i in chemicals]) / H_magnitude
        hbot = np.array([i.H(phase_bot, T, P) for i in chemicals]) / H_magnitude
        Htop = (htop * mol_top).sum()
        Hbot = (hbot * mol_bot).sum()
        JD = JacobianData(
            dEdx=dEdx,
            dHdFtop=htop,
            dHdFbot=hbot,
            dHdTtop=Ctop,
            dHdTbot=Cbot,
            split_top=self.top_split,
            split_bot=self.bottom_split,
            variable=self.specified_variable,
            value=getattr(self, self.specified_variable),
        )
        SD = StageData(
            T=T,
            top=InterstageData(Htop, mol_top, self.top_split),
            bottom=InterstageData(Hbot, mol_bot, self.bottom_split),
            feed=FeedData(H_feed / H_magnitude, mol_feed),
        )
        return JD, SD
    
    @property
    def _equilibrium_residuals_vectorized(self):
        try:
            return self._vectorized_equilibrium_residuals
        except:
            try:
                K_model = self._K_model
            except:
                self._K_model = K_model = tmo.equilibrium.PartitionCoefficients(''.join(self.phases), self.chemicals[self.partition.IDs], self.thermo)
            P = self.P
            def residuals(x):
                N_chemicals = (x.size - 1) // 2
                mol_top = x[:N_chemicals]
                mol_bottom = x[-N_chemicals:]
                T = x[N_chemicals]
                bulk_top = mol_top.sum()
                bulk_bottom = mol_bottom.sum()
                if bulk_top and bulk_bottom:
                    y = mol_top / bulk_top
                    x = mol_bottom / bulk_bottom
                    K = K_model(y, x, T, P)
                    return K * mol_bottom * bulk_top / bulk_bottom - mol_top
                elif self.B == 0 or bulk_top:
                    return -0.1 * mol_top
                else:
                    return 0.1 * mol_bottom
            
            self._vectorized_equilibrium_residuals = lambda x: np.apply_along_axis(residuals, axis=0, arr=x)
            return self._vectorized_equilibrium_residuals
    
    def _equilibrium_residuals(self, center):
        T = center.T
        mol_top = center.top.mol
        mol_bottom = center.bottom.mol
        bulk_top = mol_top.sum()
        bulk_bottom = mol_bottom.sum()
        if bulk_top and bulk_bottom:
            y = mol_top / bulk_top
            x = mol_bottom / bulk_bottom
            try:
                K_model = self._K_model
            except:
                self._K_model = K_model = tmo.equilibrium.PartitionCoefficients(''.join(self.phases), self.chemicals[self.partition.IDs], self.thermo)
            K = K_model(y, x, T, self.P)
            return K * mol_bottom * bulk_top / bulk_bottom - mol_top
        elif self.B == 0 or bulk_top:
            return -0.1 * mol_top
        else:
            return 0.1 * mol_bottom
    
    def _material_balance_residuals(self, upper, center, lower):
        mol_out = center.top.mol + center.bottom.mol
        mol_in = center.feed.mol.copy()
        inlets = []
        if upper: inlets.append(upper.bottom)
        if lower: inlets.append(lower.top)
        for inlet in inlets:
            if inlet.split:
                mol_in += (1 - inlet.split) * inlet.mol
            else:
                mol_in += inlet.mol
        return mol_out - mol_in
    
    def _energy_balance_residual(self, upper, center, lower):
        var = self.specified_variable
        if var == 'Q':
            H_out = (center.top.H + center.bottom.H) 
            H_in = center.feed.H
            inlets = []
            if upper: inlets.append(upper.bottom)
            if lower: inlets.append(lower.top)
            for i in inlets:
                H_in += (1 - i.split) * i.H
            return H_out - H_in - self.Q
        elif var == 'T':
            return self.T - center.T
        elif var == 'B':
            return center.top.mol.sum() - center.bottom.mol.sum() * self.B
        elif var == 'F':
            F = 0
            if upper:
                split = center.top.split
                if split: F += split * center.top.mol.sum()
            else:
                F += center.top.mol.sum()
            if lower:
                split = center.bottom.split
                if split: F += split * center.bottom.mol.sum()
            else:
                F += center.bottom.mol.sum()
            return self.F - F
        else:
            raise RuntimeError('unknown specification')
    
    def _get_point(self, x=None):
        IDs = self.partition.IDs
        N_chemicals = len(IDs)
        if x is None:
            x = np.zeros(N_chemicals * 2 + 1)
        top, bottom = self.partition.outs
        x[:N_chemicals] = top.imol[IDs]
        x[N_chemicals] = self.T
        x[-N_chemicals:] = bottom.imol[IDs]
        return x
    
    def _set_point(self, x):
        IDs = self.partition.IDs
        N_chemicals = len(IDs)
        top, bottom = self.partition.outs
        top.imol[IDs] = mol_top = x[:N_chemicals]
        if self.specified_variable != 'T':
            self.T = x[N_chemicals]
        bottom.imol[IDs] = mol_bot = x[-N_chemicals:]
        if self.specified_variable != 'B':
            bulk_bot = mol_bot.sum()
            if bulk_bot:
                self.B = mol_top.sum() / bulk_bot
            else:
                self.B = float('inf')
        mol_bot[mol_bot == 0] = 1-16
        self.S[:] = mol_top / mol_bot
        if self.B != 0: self.K = self.S / self.B
        for i in self.splitters: i._run()

    # %% Phenomena-based simulation
    
    @property
    def composition_sensitive(self):
        return self.phases == ('L', 'l')
    
    def _get_energy_departure_coefficient(self, stream):
        energy_variable = self._energy_variable
        if energy_variable is None: return None
        if energy_variable == 'B':
            if stream.phase != 'g': return None
            vapor, liquid = self.partition.outs
            if vapor.isempty():
                hV = self.mixture(liquid, 'h', phase='g')
            else:
                hV = vapor.h
            dHdB = hV * liquid.F_mol
            if self.top_split:
                if stream.imol is self.top_side_draw.imol:
                    split = self.top_split
                    return (self, -dHdB * split)
                elif stream.imol is self.outs[0].imol:
                    split = self.top_split
                    return (self, -dHdB * (1 - split))
                else:
                    raise ValueError('stream must be an outlet')
            elif stream.imol is self.outs[0].imol:
                return (self, -dHdB)
            else:
                raise ValueError('stream must be an outlet')
        else:
            return (self, -stream.C)
    
    def _create_energy_departure_equations(self):
        energy_variable = self._energy_variable
        if energy_variable is None: return []
        if energy_variable == 'B':
            vapor, liquid = self.partition.outs
            if vapor.isempty():
                if liquid.isempty(): 
                    raise RuntimeError('empty stage or tray')
                hV = self.mixture('h', liquid, phase='g')
            else:
                hV = vapor.h
            coeff = {self: hV * liquid.F_mol}
        else:
            coeff = {self: sum([i.C for i in self.partition.outs])}
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        if self.reaction:
            dH = self.Q - self.Hnet
        else:
            dH = self.Q + self.H_in - self.H_out
            
        return [(coeff, dH)]
    
    def _create_material_balance_equations(self, composition_sensitive):
        partition = self.partition
        chemicals = self.chemicals
        pIDs = partition.IDs
        IDs = chemicals.IDs
        self._update_separation_factors()
        if pIDs != IDs and pIDs:
            partition.IDs = IDs
            S = np.ones(chemicals.size)
            index = [IDs.index(i) for i in pIDs]
            for i, j in zip(index, self.S): S[i] = j
            pIDs = set(pIDs)
            data = self.partition.partition_data
            if data:
                top = data.get('extract_chemicals') or data.get('top_chemicals', ())
                bottom = data.get('raffinate_chemicals') or data.get('bottom_chemicals', ())
                for i in top: S[chemicals.index(i)] = inf
                for i in bottom: S[chemicals.index(i)] = 0
                pIDs.update(top)
                pIDs.update(bottom)
            for index, ID in enumerate(IDs):
                if ID in pIDs: continue
                top = partition.outs[0].mol[index]
                bottom = partition.outs[1].mol[index]
                if top:
                    if bottom:
                        S[index] =  top / bottom
                    else:
                        S[index] =  inf
                else:
                    S[index] =  0
        else:
            S = self.S.copy()
        top_split = self.top_split
        bottom_split = self.bottom_split
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        top, bottom, *_ = self.outs
        top_side_draw = self.top_side_draw
        bottom_side_draw = self.bottom_side_draw
        N = self.chemicals.size
        ones = np.ones(N)
        minus_ones = -ones
        zeros = np.zeros(N)
        
        # # Overall flows
        eq_overall = {}
        reaction = self.reaction
        if reaction: # Reactive liquid
            predetermined_flow = SparseVector.from_dict(sum_sparse_vectors([i.mol for i in fresh_inlets]), size=N)
            rhs = predetermined_flow + self.partition.dmol
            for i in self.outs: eq_overall[i] = ones
            for i in process_inlets: eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, rhs)
            )
        else:
            for i in self.outs: eq_overall[i] = ones
            for i in process_inlets:
                if i in eq_overall: del eq_overall[i]
                else: eq_overall[i] = minus_ones
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
        
        # Top to bottom flows
        eq_outs = {}
        infmask = ~np.isfinite(S)
        S[infmask] = 1
        if top_split == 1:
            if bottom_split == 1:
                eq_outs[top_side_draw] = -ones
                eq_outs[bottom_side_draw] = S
            else:
                eq_outs[top_side_draw] = -ones * (1 - bottom_split)
                eq_outs[bottom] = S
        elif bottom_split == 1:
            eq_outs[top] = coef = -ones
            eq_outs[bottom_side_draw] = S * (1 - top_split) 
            coef[infmask] = 0
        else:
            eq_outs[top] = coef = -ones * (1 - bottom_split)
            eq_outs[bottom] = S * (1 - top_split) 
            coef[infmask] = 0
        equations.append(
            (eq_outs, zeros)
        )
        # Top split flows
        if top_side_draw:
            if top_split == 1:
                eq_top_split = {
                    top: ones,
                }
            else:
                eq_top_split = {
                    top_side_draw: ones,
                    top: -top_split / (1 - top_split),
                }
            equations.append(
                (eq_top_split, zeros)
            )
        # Bottom split flows
        if bottom_side_draw:
            if bottom_split == 1:
                eq_bottom_split = {
                    bottom: ones,
                }
            else:
                eq_bottom_split = {
                    bottom_side_draw: ones,
                    bottom: -bottom_split / (1 - bottom_split),
                }
            equations.append(
                (eq_bottom_split, zeros)
            )
        return equations
    
    def _update_energy_variable(self, departure):
        phases = self.phases
        energy_variable = self._energy_variable
        if energy_variable == 'B':
            partition = self.partition
            top, bottom = partition.outs
            IDs = partition.IDs
            B = top.imol[IDs].sum() / bottom.imol[IDs].sum()
            self.B = B + departure
        elif phases == ('L', 'l'):
            self.T = T = self.T + departure
            for i in self.outs: i.T = T
        else:
            raise RuntimeError('invalid phases')
            
    def _update_composition_parameters(self):
        partition = self.partition
        data = partition.partition_data
        if data and 'K' in data: return
        partition._run_decoupled_Kgamma()
    
    def _update_net_flow_parameters(self):
        self.partition._run_decoupled_B()
    
    def _update_nonlinearities(self):
        self._update_equilibrium_variables()
        self._update_reaction_conversion()
    
    def _update_equilibrium_variables(self):
        phases = self.phases
        if phases == ('g', 'l'):
            partition = self.partition
            partition._run_decoupled_KTvle()
            T = partition.T
            for i in (*partition.outs, *self.outs): i.T = T
        elif phases == ('L', 'l'):
            pass
            # self.partition._run_lle(single_loop=True)
        else:
            raise NotImplementedError(f'K for phases {phases} is not yet implemented')
        
    def _update_reaction_conversion(self):    
        if self.reaction and self.phases == ('g', 'l'):
            self.partition._run_decoupled_reaction()
    
    # %% Graph representations and profiling
    
    @property
    def equation_node_names(self): 
        material_balances = (
            'overall_material_balance_node', 
            'separation_material_balance_node',
        )
        if self.phases == ('g', 'l'):
            phenomenode = 'vle_phenomenode'
        else: # Assume LLE
            phenomenode = 'lle_phenomenode'
        if self._energy_variable is None:
            return (
                *material_balances,
                phenomenode,
            )
        else:
            return (
                *material_balances,
                'energy_balance_node',
                phenomenode,
            )
    
    @property
    def phenomenode(self):
        return self.vle_phenomenode if self.phases == ('g', 'l') else self.lle_phenomenode
    
    def initialize_overall_material_balance_node(self):
        self.overall_material_balance_node.set_equations(
            inputs=[j for i in self.ins if (j:=i.F_node)],
            outputs=[i.F_node for i in self.outs],
        )
    
    def initialize_separation_material_balance_node(self):
        self.separation_material_balance_node.set_equations(
            outputs=[i.F_node for i in self.outs],
            inputs=[self.K_node, self.Phi_node],
        )
        
    def initialize_lle_phenomenode(self):
        intermediates = [
            i.F_node for i in self.outs 
            if hasattr(i.sink, 'lle_phenomenode')
        ]
        self.lle_phenomenode.set_equations(
            inputs=[self.T_node, *[i.F_node for i in self.ins]],
            outputs=[self.K_node, self.Phi_node, *intermediates],
            tracked_outputs=[self.K_node, self.Phi_node],
        )
    
    def initialize_vle_phenomenode(self):
        if self.specified_variable == 'T':
            self.vle_phenomenode.set_equations(
                inputs=[self.T_node, *[i.F_node for i in self.outs]],
                outputs=[self.K_node, self.Phi_node],
            )
        else:
            self.vle_phenomenode.set_equations(
                inputs=[i.F_node for i in self.outs if i.phase == 'l'],
                outputs=[self.T_node, self.K_node],
            )
            
    def initialize_energy_balance_node(self):
        self.energy_balance_node.set_equations(
            inputs=(
                self.T_node, 
                *[i.T_node for i in (*self.ins, *self.outs)],
                *[i.F_node for i in (*self.ins, *self.outs)],
                *[j for i in self.ins if (j:=i.E_node)]
            ),
            outputs=[j for i in self.outs if (j:=i.E_node)],
        )
        
    @property
    def K_node(self):
        if hasattr(self, '_K_node'): return self._K_node
        partition_data = self.partition.partition_data
        if (self.specified_variable == 'B') and (self.B == 0 or self.B == np.inf):
            var = None 
        if  (partition_data and 'K' in partition_data):
            var = None
        else:
            var = VariableNode(f"{self.node_tag}.K", lambda: self.K)
        self._K_node = var
        return var
    
    @property
    def T_node(self):
        if hasattr(self, '_T_node'): return self._T_node
        if self.specified_variable == 'T': 
            var = None
        else:
            var = VariableNode(f"{self.node_tag}.T", lambda: self.T)
        self._T_node = var
        return var 
    
    @property
    def Phi_node(self):
        if hasattr(self, '_Phi_node'): return self._Phi_node
        if self.phases == ('g', 'l'):
            if self.specified_variable == 'B':
                self._Phi_node = var = None
            else:
                self._Phi_node = var = VariableNode(f"{self.node_tag}.Phi", lambda: self.B)
        else:
            self._Phi_node = var = VariableNode(f"{self.node_tag}.Phi", lambda: self.B)
        return var
    
    def get_E_node(self, stream):
        if self.phases == ('g', 'l') and stream.phase != 'g':
            return None
        else:
            return self.E_node
    
    @property
    def E_node(self):
        if hasattr(self, '_E_node'): return self._E_node
        if self._energy_variable is None:
            var = None
        elif self.phases == ('g', 'l'):
            var = self.Phi_node
        else:
            var = self.T_node
        self._E_node = var
        return var

    
    
# %%

class PhasePartition(Unit):
    _N_ins = 1
    _N_outs = 2
    strict_infeasibility_check = False
    dmol_relaxation_factor = 0
    S_relaxation_factor = 0
    B_relaxation_factor = 0
    K_relaxation_factor = 0
    T_relaxation_factor = 0
    F_relaxation_factor = 0
    gamma_y_relaxation_factor = 0
    fgas_relaxation_factor = 0
    
    def _init(self, phases, partition_data, top_chemical=None, reaction=None):
        self.partition_data = partition_data
        self.phases = phases
        self.top_chemical = top_chemical
        self.reaction = reaction
        self.gamma_y = None
        self.fgas = None
        self.IDs = None
        self.K = None
        self.B = None
        self.T = None
        self.P = None
        self.Q = None
        self.specified_variable = None
        self.B_fallback = 1
        self.dmol = SparseVector.from_size(self.chemicals.size)
        for i, j in zip(self.outs, self.phases): i.phase = j 
        
    @property
    def x(self):
        try:
            IDs = self.IDs
            x = self.outs[1].imol[IDs]
            xsum = x.sum()
            if xsum:
                return x / xsum
            else:
                return x
        except:
            return None
        
    @property
    def y(self):
        try:
            return self.x * self.K
        except:
            return None
        
    def _get_mixture(self, linked=True):
        if linked:
            try:
                ms = self._linked_multistream 
            except:
                outs = self.outs
                for i, j in zip(self.outs, self.phases): i.phase = j 
                self._linked_multistream = ms = tmo.MultiStream.from_streams(outs)
            if self.specified_variable == 'T': ms.T = self.T
            return ms
        else:
            try:
                ms = self._unlinked_multistream
                ms.copy_like(self._get_mixture())
            except:
                self._unlinked_multistream = ms = self._get_mixture().copy()
            if self.specified_variable == 'T': ms.T = self.T
            return ms
    
    def _get_arrays(self):
        if self.gamma_y is None:
            return {'K': self.K}
        else:
            return {'K': self.K, 'gamma_y': self.gamma_y}
    
    def _set_arrays(self, IDs, **kwargs):
        IDs_last = self.IDs
        IDs = tuple(IDs)
        if IDs_last and IDs_last != IDs:
            if len(IDs_last) > len(IDs):
                index = [IDs_last.index(i) for i in IDs]
                for name, array in kwargs.items():
                    last = getattr(self, name)
                    f = getattr(PhasePartition, name + '_relaxation_factor')
                    g = 1. - f 
                    for i, j in enumerate(index):
                        last[j] = g * array[i] + f * last[j]
            else:
                self.IDs = IDs
                index = [IDs.index(i) for i in IDs_last]
                for name, array in kwargs.items():
                    last = getattr(self, name)
                    new = array.copy()
                    setattr(self, name, new)
                    for i, j in enumerate(index): new[i] = last[j]
                    f = getattr(PhasePartition, name + '_relaxation_factor')
                    g = 1. - f 
                    for i, j in enumerate(index):
                        new[j] = g * array[i] + f * new[j]
        else:
            for i, j in kwargs.items(): setattr(self, i, j)
            self.IDs = IDs
    
    def _get_activity_model(self):
        chemicals = self.chemicals
        index = chemicals.get_lle_indices(sum([i.mol for i in self.ins]).nonzero_keys())
        chemicals = chemicals.tuple
        lle_chemicals = [chemicals[i] for i in index]
        return self.thermo.Gamma(lle_chemicals), [i.ID for i in lle_chemicals], index
    
    def _get_fugacity_models(self):
        chemicals = self.chemicals
        index = chemicals.get_vle_indices(sum([i.mol for i in self.ins]).nonzero_keys())
        chemicals = chemicals.tuple
        vle_chemicals = [chemicals[i] for i in index]
        return (equilibrium.GasFugacities(vle_chemicals, thermo=self.thermo),
                equilibrium.LiquidFugacities(vle_chemicals, thermo=self.thermo), 
                [i.ID for i in vle_chemicals],
                index)
    
    def _update_fgas(self, P=None):
        F_gas, F_liq, IDs, index = self._get_fugacity_models()
        top, bottom = self.outs
        T = self.T
        y = top.mol[index]
        y_sum = y.sum()
        if y_sum: 
            y /= y_sum
        else:
            y = np.ones(y.size) / y.size
        if P is None: P = self.P
        self.fgas = F_gas.unweighted(y, T, P)
    
    def _run_decoupled_Kfgas(self, P=None):
        top, bottom = self.outs
        F_gas, F_liq, IDs, index = self._get_fugacity_models()
        if P is None: P = self.P
        T = self.T
        x = bottom.mol[index]
        x_sum = x.sum()
        if x_sum:
            x /= x_sum
        else:
            x = np.ones(x.size) / x.size
        
        try:
            fgas = self.fgas
            init = fgas is None or fgas.size != len(index)
        except:
            init = True
        if init:
            y = top.mol[index]
            y_sum = y.sum()
            if y_sum: 
                y /= y_sum
            else:
                y = np.ones(y.size) / y.size
            self.fgas = fgas = F_gas.unweighted(y, T, P)
        
        fliq = F_liq.unweighted(x, T, P)
        K = fliq / fgas 
        y = K * x
        y /= y.sum()
        fgas = F_gas.unweighted(y, T, P)
        K = fliq / fgas
        good = (x != 0) | (y != 0)
        if not good.all():
            index, = np.where(good)
            IDs = [IDs[i] for i in index]
            fgas = [fgas[i] for i in index]
            K = [K[i] for i in index]
        self._set_arrays(IDs, fgas=fgas, K=K)
    
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
            return
        if phi <= 0 or phi >= 1: return
        self.B = phi / (1 - phi)
        # TODO: set S using relaxation factor and use different separation factors for lle and vle
    
    def _run_decoupled_KTvle(self, P=None, 
                             T_relaxation_factor=None,
                             K_relaxation_factor=None): # Bubble point
        top, bottom = self.outs
        if P is not None: top.P = bottom.P = P
        if self.specified_variable == 'T':
            self._run_vle(update=False)
            for i in self.outs: i.T = self.T
        else:
            if bottom.isempty():
                if top.isempty(): return
                p = top.dew_point_at_P(P)
            else:
                p = bottom.bubble_point_at_P(P)
            # TODO: Note that solution decomposition method is bubble point
            x = p.x
            x[x == 0] = 1.
            K_new = p.y / p.x
            IDs = p.IDs
            top.imol[IDs] = p.y * top.imol[IDs].sum()
            f = self.T_relaxation_factor if T_relaxation_factor is None else T_relaxation_factor
            if self.T:
                self.T = f * self.T + (1 - f) * p.T
            else:
                self.T = p.T
            self._equilibrium_point = p
            self._set_arrays(IDs, K=K_new)
    
    def _run_decoupled_reaction(self, P=None, relaxation_factor=None):
        top, bottom = self.outs
        f = self.dmol_relaxation_factor if relaxation_factor is None else relaxation_factor
        old = self.dmol
        new = self.reaction.conversion(bottom)
        self.dmol = f * old + (1 - f) * new
    
    def _run_lle(self, P=None, update=True, top_chemical=None, single_loop=False):
        if top_chemical is None: top_chemical = self.top_chemical
        else: self.top_chemical = top_chemical
        ms = self._get_mixture(update)
        eq = ms.lle
        data = self.partition_data
        if data and 'K' in data:
            ms.phases = self.phases
            top, bottom = ms
            IDs = data['IDs']
            K = data['K']
            phi = sep.partition(
                ms, top, bottom, IDs, K, 0.5, 
                data.get('extract_chemicals') or data.get('top_chemicals'),
                data.get('raffinate_chemicals') or data.get('bottom_chemicals'),
                self.strict_infeasibility_check, 1
            )
            if phi == 1:
                self.B = np.inf
            else:
                self.B = phi / (1 - phi)
            self.K = K
            self.T = ms.T
            self.IDs = IDs
        else:
            if update:
                eq(T=ms.T, P=P or self.P, top_chemical=top_chemical, update=update, single_loop=single_loop)
                lle_chemicals, K_new, gamma_y, phi = eq._lle_chemicals, eq._K, eq._gamma_y, eq._phi
            else:
                lle_chemicals, K_new, gamma_y, phi = eq(T=ms.T, P=P or self.P, top_chemical=top_chemical, update=update)
            if phi == 1 or phi is None:
                self.B = np.inf
                self.T = ms.T
                return
            else:
                self.B = phi / (1 - phi)
            self.T = ms.T
            IDs = tuple([i.ID for i in lle_chemicals])
            self._set_arrays(IDs, K=K_new, gamma_y=gamma_y)
    
    def _run_vle(self, P=None, update=True):
        ms = self._get_mixture(update)
        kwargs = {
            self.specified_variable: getattr(self, self.specified_variable),
            'P': self.P if P is None else P,
        }
        if self.reaction: kwargs['liquid_conversion'] = self.reaction.conversion_handle(self.outs[1])
        ms.vle(**kwargs)
        index = ms.vle._index
        if self.reaction: self.dmol = ms.mol - self.feed.mol
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
            K_new = y_mol / x_mol
        else:
            K_new = np.ones(len(index)) * 1e-16
        if 'B' not in kwargs:
            if not L_total:
                self.B = inf
            else:
                self.B = V_total / L_total
        self.T = ms.T
        self._set_arrays(IDs, K=K_new)
        # TODO: Add option to set S and T using relaxation factor
    
    def _simulation_error(self):
        cache = self.T, self.B, copy(self.K), copy(self.dmol), copy(self.S), copy(self.gamma_y)
        um = getattr(self, '_unlinked_multistream', None)
        m = getattr(self, '_linked_multistream', None)
        if um is not None: self._unlinked_multistream = copy(um)
        if m is not None: self._linked_multistream = copy(m)
        error = super()._simulation_error()
        self.T, self.B, self.K, self.dmol, self.S, self.gamma_y = cache
        if um is not None: self._unlinked_multistream = um
        if m is not None: self._linked_multistream = um
        return error
    
    def _run(self):
        mixture = self._get_mixture()
        mixture.copy_like(self.feed)
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
    ...     bottom_side_draws={0: 0.673 / (1 + 0.673)},
    ...     max_attempts=10,
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
    inside_maxiter = 100
    default_max_attempts = 2
    default_maxiter = 30
    default_fallback = 'simultaneous correction', 'trust-region', 50, 1
    default_S_tolerance = 1e-6
    default_relative_S_tolerance = 1e-6
    default_algorithm = 'phenomena'
    detault_inside_out = False
    default_inner_loop_algorthm = None
    decomposition_algorithms = {
        'phenomena', 'sequential modular',
    }
    available_algorithms = {
        *decomposition_algorithms, 
        'simultaneous correction',
    }
    default_methods = {
        'phenomena': 'wegstein',
        'simultaneous correction': 'trust-region',
    }
    method_options = {
        'fixed-point': {},
        'wegstein': {'lb': 1, 'ub': 4, 'exp': 0.5}
    }
    auxiliary_unit_names = (
        'stages',
    )
    _side_draw_names = ('top_side_draws', 'bottom_side_draws')
    
    # Inside-out surrogate model
    TSurrogate = RBFInterpolator
    KSurrogate = Akima1DInterpolator 
    hSurrogate = Akima1DInterpolator 
    T_surrogate_options = {}
    K_surrogate_options = {'method': 'makima', 'extrapolate': True}
    h_surrogate_options = {'method': 'makima', 'extrapolate': True}
    
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
    
    def __init__(self,  ID='', ins=None, outs=(), thermo=None, stages=None, **kwargs):
        if stages is None:
            if 'feed_stages' in kwargs: self._N_ins = len(kwargs['feed_stages'])
            top_side_draws, bottom_side_draws = self._side_draw_names
            N_outs = 2
            if top_side_draws in kwargs: N_outs += len(kwargs[top_side_draws]) 
            if bottom_side_draws in kwargs: N_outs += len(kwargs[bottom_side_draws]) 
            self._N_outs = N_outs
            Unit.__init__(self, ID, ins, outs, thermo, **kwargs)
        else:
            ins = []
            outs = []
            top_side_draws_outs = []
            bottom_side_draws_outs = []
            stages_set = set(stages)
            top_side_draws = {}
            bottom_side_draws = {}
            feed_stages = []
            first_stage = stages[0]
            phases = first_stage.phases
            stage_specifications = {}
            stage_reactions = {}
            self._load_thermo(thermo or first_stage.thermo)
            for n, stage in enumerate(stages):
                for s in stage.ins:
                    if s.source not in stages_set: 
                        sp = s.proxy()
                        sp._source = s._source
                        ins.append(sp)
                        feed_stages.append(n)
                top, bottom, *other = stage.outs
                if stage.top_split:
                    s = other[0]
                    sp = s.proxy()
                    sp._sink = s._sink
                    top_side_draws_outs.append(sp)
                    top_side_draws[n] = stage.top_split
                if stage.bottom_split:
                    s = other[-1]
                    sp = s.proxy()
                    sp._sink = s._sink
                    bottom_side_draws_outs.append(sp)
                    bottom_side_draws[n] = stage.bottom_split
                if top.sink not in stages_set: 
                    sp = top.proxy()
                    sp._sink = top._sink
                    outs.append(sp)
                if bottom.sink not in stages_set: 
                    sp = bottom.proxy()
                    sp._sink = bottom._sink
                    outs.append(sp)
                specified_variable = stage.specified_variable
                if specified_variable == 'B': 
                    stage_specifications[n] = ('Boilup', stage.B)
                elif specified_variable == 'T': 
                    stage_specifications[n] = ('Temperature', stage.T)
                elif specified_variable == 'Q':
                    stage_specifications[n] = ('Duty', stage.Q)
                elif specified_variable == 'F':
                    stage_specifications[n] = ('Flow', stage.F)
                if stage.reaction is not None:
                    stage_reactions[n] = stage.reaction
            outs = [*outs, *top_side_draws_outs, *bottom_side_draws_outs]
            self._N_ins = len(ins)
            self._N_outs = len(outs)
            Unit.__init__(self, ID, ins, outs, thermo, 
                stage_specifications=stage_specifications,
                stage_reactions=stage_reactions,
                feed_stages=feed_stages,
                bottom_side_draws=bottom_side_draws,
                top_side_draws=top_side_draws,
                stages=stages,
                phases=phases,
                **kwargs
            )
    
    def _init(self,
            N_stages=None, 
            stages=None,
            top_side_draws=None,
            bottom_side_draws=None, 
            feed_stages=None, 
            phases=None, 
            P=101325, 
            T=None,
            stage_specifications=None, 
            stage_reactions=None,
            partition_data=None, 
            top_chemical=None, 
            use_cache=None,
            collapsed_init=False,
            algorithm=None,
            method=None,
            maxiter=None,
            max_attempts=None,
            inside_out=None,
            vle_decomposition=None,
        ):
        # For VLE look for best published algorithm (don't try simple methods that fail often)
        if N_stages is None: N_stages = len(stages)
        if phases is None: phases = ('g', 'l')
        if feed_stages is None: feed_stages = (0, -1)
        if stage_specifications is None: stage_specifications = {}
        elif not isinstance(stage_specifications, dict): stage_specifications = dict(stage_specifications)
        if T is not None: 
            for i in range(N_stages):
                if i in stage_specifications: continue
                stage_specifications[i] = ('Temperature', T)
        if stage_reactions is None: stage_reactions = {}
        elif not isinstance(stage_reactions, dict): stage_reactions = dict(stage_reactions)
        if top_side_draws is None: top_side_draws = {}
        elif not isinstance(top_side_draws, dict): top_side_draws = dict(top_side_draws)
        if bottom_side_draws is None: bottom_side_draws = {}
        elif not isinstance(bottom_side_draws, dict): bottom_side_draws = dict(bottom_side_draws)
        if partition_data is None: partition_data = {}
        self.multi_stream = tmo.MultiStream(None, P=P, phases=phases, thermo=self.thermo)
        self.N_stages = N_stages
        self.P = P
        self.T = T
        self.phases = phases = self.multi_stream.phases # Corrected order
        self._has_vle = 'g' in phases
        self._has_lle = 'L' in phases
        self._top_split = top_splits = np.zeros(N_stages)
        self._bottom_split = bottom_splits = np.zeros(N_stages)
        if stages is None:
            top_mark = 2 + len(top_side_draws)
            tsd_iter = iter(self.outs[2:top_mark])
            bsd_iter = iter(self.outs[top_mark:])
            last_stage = None
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
                    try:
                        outs.append(next(bsd_iter))
                    except:
                        breakpoint()
                    bottom_split = bottom_side_draws[i]
                    bottom_splits[i] = bottom_split
                else: 
                    bottom_split = 0
                try:
                    new_stage = self.auxiliary(
                        'stages', StageEquilibrium, phases=phases,
                        ins=feed,
                        outs=outs,
                        partition_data=partition_data,
                        top_split=top_split,
                        bottom_split=bottom_split,
                    )
                except:
                    breakpoint()
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
            #: dict[int, tuple(str, float)] Specifications for VLE by stage
            self.stage_specifications = stage_specifications
            for i, (name, value) in stage_specifications.items():
                stages[i].specify_variables(**{name: value})
            self.stage_reactions = stage_reactions
            for i, reaction in stage_reactions.items():
                stages[i].reaction = reaction
        else:
            self.stage_specifications = stage_specifications
            self.stage_reactions = stage_reactions
            self.stages = stages
            top_splits = np.zeros(N_stages)
            bottom_splits = top_splits.copy()
            for i, j in top_side_draws.items(): top_splits[i] = j
            for i, j in bottom_side_draws.items(): bottom_splits[i] = j
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
            
        #: [int] Maximum number of iterations.
        self.maxiter = self.default_maxiter if maxiter is None else maxiter
        
        #: [int] Maximum number of attempts.
        self.max_attempts = self.default_max_attempts if max_attempts is None else max_attempts
        
        #: tuple[str, str] Fallback algorithm and method.
        self.fallback = self.default_fallback

        #: [float] Separation factor tolerance
        self.S_tolerance = self.default_S_tolerance

        #: [float] Relative separation factor tolerance
        self.relative_S_tolerance = self.default_relative_S_tolerance
        
        self.use_cache = True if use_cache else False
        
        self.collapsed_init = collapsed_init
        
        self.algorithm = self.default_algorithm if algorithm is None else algorithm
        
        self.method = self.default_methods[self.algorithm] if method is None else method
        
        self.inside_out = self.detault_inside_out if inside_out is None else inside_out
        
        self.vle_decomposition = vle_decomposition
    
    # %% Optimization
    
    def _get_point(self):
        N_chemicals = self._N_chemicals
        N_stages = self.N_stages
        N_variables = 2 * N_chemicals + 1
        x = np.zeros([N_stages, N_variables])
        for i, stage in enumerate(self.stages): 
            stage._get_point(x[i])
        return x
    
    def _set_point(self, x):
        for i, stage in enumerate(self.stages): 
            stage._set_point(x[i])
    
    def _objective(self, x):
        mask = x < 1e-12
        if mask.any():
            lb = x[mask].min()
            ub = x[mask].max()
            distance = ub - lb
            if distance:
                values = x[mask]
                x[mask] = (values - lb) * 1e-12 / distance + 1e-16
            else:
                x[mask] = 1e-12
        return self._net_residual(self._residuals(x))
    
    def _net_residual(self, residuals):
        net_residual = (residuals * residuals).sum()
        self._last_residuals = (residuals, net_residual)
        return net_residual
    
    def _residuals(self, x):
        N_chemicals = self._N_chemicals
        N_stages, N_variables = x.shape
        stages = self.stages
        residuals = np.zeros(N_variables * N_stages) # H, Mi, Ei
        H_magnitude = self._H_magnitude
        stage_data = [
            stage._stage_data(xi, H_feed, mol_feed, H_magnitude)
            for xi, stage, H_feed, mol_feed
            in zip(x, stages, self.feed_enthalpies, self.feed_flows)
        ]
        center = stage_data[0]
        lower = stage_data[1]
        stage = stages[0]
        H_index = 0
        M_slice = slice(1, N_chemicals + 1)
        E_slice = slice(N_chemicals + 1, None)
        residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
        residuals[0, H_index] = stage._energy_balance_residual(None, center, lower)
        residuals[0, M_slice] = stage._material_balance_residuals(None, center, lower)
        stage = stages[1]
        for i in range(2, N_stages): 
            upper = center
            center = lower
            lower = stage_data[i]
            ilast = i-1
            residuals[ilast, H_index] = stage._energy_balance_residual(upper, center, lower)
            residuals[ilast, M_slice] = stage._material_balance_residuals(upper, center, lower)
            residuals[ilast, E_slice] = stage._equilibrium_residuals(center)
            stage = stages[i]
        upper = center
        center = lower
        residuals[i, H_index] = stage._energy_balance_residual(upper, center, None)
        residuals[i, M_slice] = stage._material_balance_residuals(upper, center, None)
        residuals[i, E_slice] = stage._equilibrium_residuals(center)
        return residuals
    
    def _jacobian(self, x): # returns diagonal blocks
        N_chemicals = self._N_chemicals
        N_stages, N_variables = x.shape
        stages = self.stages
        H_magnitude = self._H_magnitude
        jacobian_data = [stage._jacobian_data(xi, H_magnitude) for xi, stage in zip(x, stages)]
        return jacobian_blocks(jacobian_data, N_stages, N_chemicals, N_variables)
        
    def _correction(self, x):
        N_chemicals = self._N_chemicals
        N_stages, N_variables = x.shape
        jacobian_data = []
        stage_data = []       
        stages = self.stages
        H_magnitude = self._H_magnitude
        for i, (xi, stage, H_feed, mol_feed) in enumerate(zip(x, stages, self.feed_enthalpies, self.feed_flows)):
            JD, SD = stage._simultaneous_correction_data(xi, H_feed, mol_feed, H_magnitude)
            jacobian_data.append(JD)
            stage_data.append(SD)
        center = stage_data[0]
        lower = stage_data[1]
        stage = stages[0]
        H_index = 0
        M_slice = slice(1, N_chemicals + 1)
        E_slice = slice(N_chemicals + 1, None)
        residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
        residuals[0, H_index] = stage._energy_balance_residual(None, center, lower)
        residuals[0, M_slice] = stage._material_balance_residuals(None, center, lower)
        stage = stages[1]
        for i in range(2, N_stages): 
            upper = center
            center = lower
            lower = stage_data[i]
            ilast = i-1
            residuals[ilast, H_index] = stage._energy_balance_residual(upper, center, lower)
            residuals[ilast, M_slice] = stage._material_balance_residuals(upper, center, lower)
            residuals[ilast, E_slice] = stage._equilibrium_residuals(center)
            stage = stages[i]
        upper = center
        center = lower
        residuals[i, H_index] = stage._energy_balance_residual(upper, center, None)
        residuals[i, M_slice] = stage._material_balance_residuals(upper, center, None)
        residuals[i, E_slice] = stage._equilibrium_residuals(center)
        A, B, C = jacobian_blocks(jacobian_data, N_stages, N_chemicals, N_variables)
        correction = MESH.solve_block_tridiagonal_matrix(A, B, C, residuals)
        return correction
    
    def _simultaneous_correction_iter(self, x): 
        # Simple line search, fsolve (which uses trust-region) is better
        try:
            residuals, residual = self._last_residuals
        except:
            N_chemicals = self._N_chemicals
            N_stages, N_variables = x.shape
            jacobian_data = []
            stage_data = []       
            stages = self.stages
            H_magnitude = self._H_magnitude
            residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
            for i, (xi, stage, H_feed, mol_feed) in enumerate(zip(x, stages, self.feed_enthalpies, self.feed_flows)):
                JD, SD = stage._simultaneous_correction_data(xi, H_feed, mol_feed, H_magnitude)
                jacobian_data.append(JD)
                stage_data.append(SD)
            center = stage_data[0]
            lower = stage_data[1]
            stage = stages[0]
            H_index = 0
            M_slice = slice(1, N_chemicals + 1)
            E_slice = slice(N_chemicals + 1, None)
            residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
            residuals[0, H_index] = stage._energy_balance_residual(None, center, lower)
            residuals[0, M_slice] = stage._material_balance_residuals(None, center, lower)
            stage = stages[1]
            for i in range(2, N_stages): 
                upper = center
                center = lower
                lower = stage_data[i]
                ilast = i-1
                residuals[ilast, H_index] = stage._energy_balance_residual(upper, center, lower)
                residuals[ilast, M_slice] = stage._material_balance_residuals(upper, center, lower)
                residuals[ilast, E_slice] = stage._equilibrium_residuals(center)
                stage = stages[i]
            upper = center
            center = lower
            residuals[i, H_index] = stage._energy_balance_residual(upper, center, None)
            residuals[i, M_slice] = stage._material_balance_residuals(upper, center, None)
            residuals[i, E_slice] = stage._equilibrium_residuals(center)
            A, B, C = jacobian_blocks(jacobian_data, N_stages, N_chemicals, N_variables)
            residual = self._net_residual(residuals)
        else:
            A, B, C = self._jacobian(x)
        correction = -MESH.solve_block_tridiagonal_matrix(A, B, C, residuals)
        try:
            tguess = self._tguess
        except:
            self._tguess = tguess = 1
        try:
            tguess, x, new_residual = flx.inexact_line_search(self._objective, x, correction, fx=residual, tguess=tguess)
        except:
            self._last_residuals = residuals, residual # Reset to original
            return 0, x, residual
        residuals, residual = self._last_residuals
        self._tguess = tguess
        return tguess, x, residual
    
    # %% Decoupled phenomena equation oriented simulation
    
    @property
    def composition_sensitive(self):
        return self._has_lle
    
    def _get_energy_departure_coefficient(self, stream):
        if self._has_vle:
            vapor, liquid = self.outs
            if stream.imol is vapor.imol:
                if vapor.isempty():
                    with liquid.temporary_phase('g'): coeff = liquid.H
                else:
                    coeff = -vapor.h * liquid.F_mol
        else:
            coeff = -stream.C
        return (self, coeff)
    
    def _create_energy_departure_equations(self):
        # Ll: C1dT1 - Ce2*dT2 - Cr0*dT0 - hv2*L2*dB2 = Q1 - H_out + H_in
        # gl: hV1*L1*dB1 - hv2*L2*dB2 - Ce2*dT2 - Cr0*dT0 = Q1 + H_in - H_out
        phases = self.phases
        if phases == ('g', 'l'):
            vapor, liquid = self.outs
            coeff = {}
            if vapor.isempty():
                with liquid.temporary_phase('g'): coeff[self] = liquid.H
            else:
                coeff[self] = vapor.h * liquid.F_mol
        elif phases == ('L', 'l'):
            coeff = {self: sum([i.C for i in self.outs])}
        else:
            raise RuntimeError('invalid phases')
        for i in self.ins: i._update_energy_departure_coefficient(coeff)
        if self.stage_reactions:
            return [(coeff, sum([i.Q for i in self.stages]) - self.Hnet)]
        else:
            return [(coeff, self.H_in - self.H_out + sum([(i.Hnet if i.Q is None else i.Q) for i in self.stages]))]
    
    def _create_material_balance_equations(self, composition_sensitive):
        top, bottom = self.outs
        try:
            B = self.B
            K = self.K
        except:
            if bottom.isempty():
                self.B = B = np.inf
                self.K = K = 1e16 * np.ones(self.chemicals.size)
            elif top.isempty():
                self.K = K = np.zeros(self.chemicals.size)
                self.B = B = 0
            else:
                top_mol = top.mol.to_array()
                bottom_mol = bottom.mol.to_array()
                Ftop = top_mol.sum()
                Fbot = bottom_mol.sum()
                y = top_mol / Ftop
                x = bottom_mol / Fbot
                x[x <= 0] = 1e-16
                self.K = K = y / x
                self.B = B = Ftop / Fbot
        
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        top, bottom, *_ = self.outs
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {}
        for i in self.outs: eq_overall[i] = ones
        for i in process_inlets: 
            if i in eq_overall:
                del eq_overall[i]
            else:
                eq_overall[i] = minus_ones
        if self.stage_reactions:
            partitions = self.partitions
            flows = [i.mol for i in fresh_inlets] + [partitions[i].dmol for i in self.stage_reactions]
            equations.append(
                (eq_overall, sum(flows, zeros))
            )
        else:
            equations.append(
                (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
            )
        
        # Top to bottom flows
        eq_outs = {}
        if B == np.inf:
            eq_outs[bottom] = ones
        elif B == 0:
            eq_outs[top] = ones
        else:
            eq_outs[top] = ones
            eq_outs[bottom] = -(K * B)
        equations.append(
            (eq_outs, zeros)
        )
        return equations
    
    def _update_composition_parameters(self):
        for i in self.partitions: 
            if 'K' in i.partition_data: continue
            i._run_decoupled_Kgamma()
    
    def _update_net_flow_parameters(self):
        for i in self.partitions: i._run_decoupled_B()
    
    def _update_nonlinearities(self):
        if self._has_vle:
            for i in self.stages: i._update_nonlinearities()
        elif self._has_lle:
            pass
            # self.update_pseudo_lle()
    
    def _update_energy_variable(self, departure):
        phases = self.phases
        if phases == ('g', 'l'):
            if not hasattr(self, 'B'):
                top, bottom = self.outs
                if bottom.isempty():
                    self.B = np.inf
                    self.K = 1e16 * np.ones(self.chemicals.size)
                elif top.isempty():
                    self.K = np.zeros(self.chemicals.size)
                    self.B = 0
                else:
                    top_mol = top.mol.to_array()
                    bottom_mol = bottom.mol.to_array()
                    Ftop = top_mol.sum()
                    Fbot = bottom_mol.sum()
                    y = top_mol / Ftop
                    x = bottom_mol / Fbot
                    x[x <= 0] = 1e-16
                    self.K = y / x
                    self.B = Ftop / Fbot
            self.B += departure
        elif phases == ('L', 'l'):
            for i in self.outs: i.T += departure
        else:
            raise RuntimeError('invalid phases')
    
    # %% Supportive utilities
    
    def correct_overall_mass_balance(self):
        outmol = sum([i.mol for i in self.outs])
        inmol = sum([i.mol for i in self.ins])
        stage_reactions = self.stage_reactions
        if stage_reactions:
            partitions = self.partitions
            inmol += sum([partitions[i].dmol for i in stage_reactions])
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
                sum([i.mol for i in stage.ins],
                    -sum([i.mol for i in stage.outs], -stage.partition.dmol))
            )
        return pd.DataFrame(errors, columns=IDs)
    
    # %% Simulation    
    
    def default_vle_decomposition(self):
        K = np.mean([i.K for i in self.stages], axis=0)
        mol = self.feed_flows.sum(axis=0)
        z = mol / mol.sum()
        if equilibrium.stable_phase(K, z):
            self.vle_decomposition = 'sum rates'
        else:
            self.vle_decomposition = 'bubble point'
        
    def _run(self):
        if all([i.isempty() for i in self.ins]): 
            for i in self.outs: i.empty()
            return
        try:
            separation_factors = self.hot_start()
            algorithm = self.algorithm
            method = self.method
            if algorithm in self.decomposition_algorithms:
                options = dict(
                    maxiter=self.maxiter,
                    xtol=self.S_tolerance,
                    rtol=self.relative_S_tolerance,
                    **self.method_options[method],
                )
                if method == 'fixed-point':
                    solver = flx.fixed_point
                elif method == 'wegstein':
                    solver = flx.wegstein
                else:
                    raise ValueError(f'invalid method {method!r}')
                f = self._inside_out_iter if self.inside_out else self._iter
                self.attempt = 0
                self.bottom_flows = None
                if self.vle_decomposition is None:
                    self.default_vle_decomposition()
                for n in range(self.max_attempts):
                    self.attempt = n
                    self.iter = 0
                    try:
                        separation_factors = solver(f, separation_factors, **options)
                    except:
                        for i in self.stages: i._run()
                        for i in reversed(self.stages): i._run()
                        separation_factors = np.array([i.S for i in self._S_stages])
                        if self.fallback and self.fallback[0] != self.algorithm:
                            original = self.algorithm, self.method, self.maxiter, self.max_attempts
                            self.algorithm, self.method, self.maxiter, self.max_attempts = self.fallback
                            try:
                                self._run()
                            finally:
                                self.algorithm, self.method, self.maxiter, self.max_attempts = original
                    else:
                        break
            elif algorithm == 'simultaneous correction':
                F_mol = self.feed_flows.sum()
                self._H_magnitude = sum([i.Hvap for i in self.ins]) / F_mol
                x = self._get_point()
                shape = x.shape
                x = x.flatten()
                f = lambda x: self._residuals(x.reshape(shape)).flatten()
                jac = lambda x: MESH.create_block_tridiagonal_matrix(*self._jacobian(x.reshape(shape)))
                x, *self._simultaneous_correction_info = fsolve(f, x, fprime=jac, full_output=True, maxfev=self.maxiter, xtol=self.S_tolerance * F_mol)
                self._set_point(x.reshape(shape))
            else:
                raise RuntimeError(
                    f'invalid algorithm {algorithm!r}, only {self.available_algorithms} are allowed'
                )
            # self.correct_overall_mass_balance()
        except Exception as e:
            if self.use_cache:
                self.use_cache = False
                try:
                    self._run()
                finally:
                    self.use_cache = True
            else:
                raise e
    
    def _phenomena_iter(self, separation_factors):
        if self._has_vle:
            decomp = self.vle_decomposition
            if decomp == 'bubble point':
                self.update_flow_rates(separation_factors, update_B=True)
                self.update_bubble_point()
                self.update_energy_balance_phase_ratio_departures()
                for i in self.stages: i._update_separation_factors()
            elif decomp == 'sum rates':
                self.update_pseudo_vle(separation_factors)
                self.update_energy_balance_temperatures()
            else:
                raise NotImplementedError(f'{decomp!r} decomposition not implemented')
        elif self._has_lle:
            self.update_pseudo_lle(separation_factors)
            self.update_energy_balance_temperatures()
        else:
            raise RuntimeError('unknown equilibrium phenomena')
        if not hasattr(self, 'deltas'): self.deltas = []
        x = lambda stream: np.array([*stream.mol, stream.H])
        self.deltas.append(
            [x(i.outs[0]) - sum([x(i) for i in i.ins if i.phase=='g']) for i in self.partitions]
        )
        return np.array([i.S for i in self._S_stages])
    
    def _inside_out_iter(self, separation_factors):
        self.iter += 1
        self.inside_iter = 0
        separation_factors = self._iter(separation_factors)
        if self._update_surrogate_models():
            return flx.fixed_point(
                self._surrogate_iter, 
                separation_factors, 
                maxiter=self.inside_maxiter,
                xtol=self.S_tolerance,
                rtol=self.relative_S_tolerance, 
                checkconvergence=False,
                checkiter=False,
            )
        else:
            return separation_factors
    
    def _update_surrogate_models(self):
        partitions = self.partitions
        mixture = self.mixture
        range_stages = range(self.N_stages)
        P = self.P
        Ts = np.array([i.T for i in partitions])
        for i in range(self.N_stages-1):
            if Ts[i+1] - Ts[i] < 0: return False
        xs = np.array([i.x for i in partitions])
        ys = np.array([i.y for i in partitions])
        Ks = np.array([i.K for i in partitions])
        hVs = np.array([
            mixture.H(mol=ys[i], phase='g', T=Ts[i], P=P)
            for i in range_stages
        ])
        hLs = np.array([
            mixture.H(mol=xs[i], phase='l', T=Ts[i], P=P)
            for i in range_stages
        ])
        self.Tlb = Ts.min()
        self.Tub = Ts.max()
        self.T_surrogate = self.TSurrogate(xs[:, :-1], Ts, **self.T_surrogate_options)
        self.K_surrogate = self.KSurrogate(Ts, Ks, **self.K_surrogate_options)
        self.hV_surrogate = self.hSurrogate(Ts, hVs, **self.h_surrogate_options)
        self.hL_surrogate = self.hSurrogate(Ts, hLs, **self.h_surrogate_options)
        # import matplotlib.pyplot as plt
        # xs = np.linspace(self.Tlb, self.Tub, num=100)
        # fig, ax = plt.subplots()
        # ax.plot(Ts, hVs, "o", label="data")
        # ax.plot(xs, self.hV_surrogate(xs), label='surrogate')
        # ax.legend()
        # fig.show()
        return True
    
    def _sequential_iter(self, separation_factors):
        self.update_flow_rates(separation_factors)
        for i in self.stages: i._run()
        for i in reversed(self.stages): i._run()
        return np.array([i.S for i in self._S_stages])
    
    def _phenomena_surrogate_iter(self, separation_factors):
        feed_flows, *args = self._iter_args
        stages = self.stages
        S_stages = self._S_stages
        for i, j in zip(S_stages, separation_factors): i.S = j
        separation_factors = np.array([i.S for i in stages])
        bottom_flows = MESH.bottom_flow_rates(separation_factors, feed_flows, *args)
        N_stages = self.N_stages
        range_stages = range(N_stages)
        top_flows = MESH.top_flows_mass_balance(
            bottom_flows, feed_flows, self._asplit_left, self._bsplit_left, N_stages
        )
        for i in range_stages:
            t = top_flows[i]
            mask = t < 0
            if mask.any():
                bulk_t = t.sum()
                t[mask] = 0
                dummy = t.sum()
                if dummy: t *= bulk_t / dummy
            b = bottom_flows[i]
            mask = b < 0
            if mask.any():
                bulk_b = b.sum()
                b[mask] = 0
                dummy = b.sum()
                if dummy: b *= bulk_b / b.sum()
        V = top_flows.sum(axis=1)
        L = bottom_flows.sum(axis=1, keepdims=True)
        L[L == 0] = 1
        xs = bottom_flows / L
        L = L[:, 0]
        Ts = self.T_surrogate(xs[:, :-1])
        Ts[Ts < self.Tlb] = self.Tlb
        Ts[Ts > self.Tub] = self.Tub
        for i, j in enumerate(stages):
            if j.specified_variable == 'T': Ts[i] = j.T
        Ks = self.K_surrogate(Ts)
        hV = self.hV_surrogate(Ts)
        hL = self.hL_surrogate(Ts)
        specification_index = [
            i for i, j in enumerate(self.stages)
            if j._energy_variable
        ]
        feed_enthalpies = self.feed_enthalpies
        specification_index = np.array(specification_index, dtype=int)
        dB = MESH.phase_ratio_departures(
            L, V, hL, hV, 
            self._asplit_1, 
            self._asplit_left,
            self._bsplit_left,
            N_stages,
            specification_index,
            feed_enthalpies,
        )
        Bs = np.array([i.B for i in self.stages])
        Bs_spec = Bs[specification_index]
        Bs += dB * 0.1
        Bs[specification_index] = Bs_spec
        for i, j in zip(stages, Bs): i.B = j
        mask = Bs != 0
        S = Bs[mask][:, None] * Ks[mask]
        for i, j in zip(S_stages, S): i.S = S
        print('-------')
        print(self.iter, self.inside_iter)
        print('-------')
        print(Ts)
        print(hL)
        print(hV)
        # breakpoint()
        return S
    
    def _surrogate_iter(self, separation_factors=None):
        self.inside_iter += 1
        algorithm = self.algorithm
        if algorithm == 'phenomena':
            separation_factors = self._phenomena_surrogate_iter(separation_factors)
        else:
            raise RuntimeError(f'invalid algorithm {algorithm!r} for surrogate model')
        return separation_factors
    
    def _iter(self, separation_factors=None):
        self.iter += 1
        algorithm = self.algorithm
        if algorithm == 'phenomena':
            separation_factors = self._phenomena_iter(separation_factors)
        elif algorithm == 'sequential modular':
            separation_factors = self._sequential_iter(separation_factors)
        else:
            raise RuntimeError(f'invalid algorithm {algorithm!r}')
        return separation_factors
    
    # %% Initial guess
    
    def _hot_start_phase_ratios_iter(self, 
            top_flow_rates, *args
        ):
        bottom_flow_rates = MESH.hot_start_bottom_flow_rates(
            top_flow_rates, *args
        )
        top_flow_rates = MESH.hot_start_top_flow_rates(
            bottom_flow_rates, *args
        )
        return top_flow_rates
        
    def hot_start_phase_ratios(self):
        stages = self.stages
        stage_index = []
        phase_ratios = []
        for i in list(self.stage_specifications):
            p = stages[i].partition
            if p.specified_variable != 'B': continue 
            stage_index.append(i)
            phase_ratios.append(p.B)
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
            top_flow_rates, args=args, xtol=self.relative_S_tolerance,
            checkiter=False,
        )
        bottom_flow_rates = MESH.hot_start_bottom_flow_rates(
            top_flow_rates, *args
        )
        bf = bottom_flow_rates.sum(axis=1)
        bf[bf == 0] = 1e-32
        return top_flow_rates.sum(axis=1) / bf
    
    def hot_start_collapsed_stages(self,
            all_stages, feed_stages, stage_specifications,
            top_side_draws, bottom_side_draws,
        ):
        last = 0
        for i in sorted(all_stages):
            if i == last + 1: continue
            all_stages.add(i)
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
                partition.fgas = collapsed_partition.fgas
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
        if self.stage_reactions:
            nonzero = set()
            for rxn in self.stage_reactions.values():
                nonzero.update(rxn.stoichiometry.nonzero_keys())
            all_IDs = set(IDs)
            for i in nonzero:
                ID = chemicals.IDs[i]
                if ID not in all_IDs:
                    IDs.append(ID)
        self._IDs = IDs = tuple(IDs)
        self._N_chemicals = N_chemicals = len(IDs)
        self._S_stages = [i for i in stages if i.specified_variable != 'B' or i.B != 0]
        self._NS_stages = len(self._S_stages)
        self._update_index = index = ms.chemicals.get_index(IDs)
        self.feed_flows = feed_flows = np.zeros([N_stages, N_chemicals])
        self.feed_enthalpies = feed_enthalpies = np.zeros(N_stages)
        for feed, stage in zip(feeds, feed_stages):
            feed_flows[stage, :] += feed.mol[index]
            feed_enthalpies[stage] += feed.H
        self.total_feed_flows = feed_flows.sum(axis=1)
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
                phi = data.get('phi', 0.5)
                if K.ndim == 2:
                    data['phi'] = phi = sep.partition(
                        ms, top, bottom, IDs, K.mean(axis=0), phi,
                        top_chemicals, bottom_chemicals
                    )
                    B = inf if phi == 1 else phi / (1 - phi)
                    T = ms.T
                    for i, Ki in zip(partitions, K): 
                        if i.specified_variable != 'B': i.B = B
                        i.T = T
                        i.K = Ki
                else:
                    data['phi'] = phi = sep.partition(ms, top, bottom, IDs, K, phi,
                                                      top_chemicals, bottom_chemicals)
                    B = inf if phi == 1 else phi / (1 - phi)
                    T = ms.T
                    for i in partitions: 
                        if i.specified_variable != 'B': i.B = B
                        i.T = T
                        i.K = K
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
                    i.K = K
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
                    z = bp.z
                    z[z == 0] = 1.
                    x = dp.x
                    x[x == 0] = 1.
                    K_dew = dp.z / dp.x
                    K_bubble = bp.y / bp.z
                    dK_stage = (K_bubble - K_dew) / N_stages
                    for i, B in enumerate(phase_ratios):
                        partition = partitions[i]
                        if partition.specified_variable != 'B': partition.B = B
                        partition.T = T = T_top + i * dT_stage
                        partition.K = K_dew + i * dK_stage
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
                    B = phi / (1 - phi)
                    T = ms.T
                    for partition in partitions:
                        partition.T = T
                        partition.B = B
                        for i in partition.outs: i.T = T
                        partition.K = K
                        partition.fgas = P * y_mol
                        for s in partition.outs: s.empty()
            N_chemicals = len(index)
        if top_chemicals:
            top_side_draws = self.top_side_draws
            n = len(top_chemicals)
            b = np.ones([N_stages, n])
            c = self._asplit_1[1:]
            d = np.zeros([N_stages, n])
            for feed, stage in zip(feeds, feed_stages):
                d[stage] += feed.imol[top_chemicals]
            top_flow_rates = MESH.solve_RBDMA(b, c, d)
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
            bottom_flow_rates = MESH.solve_left_bidiagonal_matrix(a, b, d)
            for partition, b in zip(partitions, bottom_flow_rates):
                partition.outs[1].imol[bottom_chemicals] = b
        if top_chemicals or bottom_chemicals:
            for i in stages:
                for s in i.splitters: s._run()
        for i in partitions: i.IDs = IDs
        self.interpolate_missing_variables()
        for i in self.stages: i._init_separation_factors()
        return np.array([i.S for i in self._S_stages])
    
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
        N_chemicals = self._N_chemicals
        for i in range(N_stages):
            partition = partitions[i]
            B = partition.B
            T = partition.T
            K = partition.K
            if B is None or K is None or K.size != N_chemicals: continue
            index.append(i)
            Bs.append(B)
            Ks.append(K)
            Ts.append(T)
            if lle: gamma_y.append(partition.gamma_y)
        N_ok = len(index)
        if len(index) != N_stages:
            if N_ok > 1:
                neighbors = MESH.get_neighbors(index, size=N_stages)
                Bs = MESH.fillmissing(neighbors, MESH.expand(Bs, index, N_stages))
                Ts = MESH.fillmissing(neighbors, MESH.expand(Ts, index, N_stages))
                N_chemicals = self._N_chemicals
                all_Ks = np.zeros([N_stages, N_chemicals])
                if lle: all_gamma_y = all_Ks.copy()
                for i in range(N_chemicals):
                    all_Ks[:, i] = MESH.fillmissing(
                        neighbors, 
                        MESH.expand([stage[i] for stage in Ks], index, N_stages)
                    )
                    if not lle: continue
                    all_gamma_y[:, i] = MESH.fillmissing(
                        neighbors, 
                        MESH.expand([stage[i] for stage in gamma_y], index, N_stages)
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
                if partition.specified_variable != 'B': partition.B = Bs[i]
                partition.K = Ks[i]
                if lle: partition.gamma_y = gamma_y[i]
    
    # %% Energy balance convergence
    
    def get_energy_balance_temperature_departures(self):
        partitions = self.partitions
        T_specified = [i.specified_variable == 'T' for i in partitions]
        N_stages = self.N_stages
        if any(T_specified):
            start = 0
            Cl = np.zeros(N_stages)
            Cv = Cl.copy()
            Hv = Cl.copy()
            Hl = Cl.copy()
            dT = Cl.copy()
            for i, p in enumerate(partitions):
                if T_specified[i]:
                    top, bottom = p.outs
                    Hl[i] = bottom.H
                    Hv[i] = top.H
                    Cl[i] = bottom.C
                    Cv[i] = top.C
                else:
                    end = i + 1
                    index = slice(start, end)
                    dT[index] = MESH.temperature_departures(
                        Cv[index], Cl[index], Hv[index], Hl[index], 
                        self._asplit_left[index], 
                        self._bsplit_left[index],
                        end - start, self.feed_enthalpies[index],
                    )
                    start = end
        else:
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
            dTs = MESH.temperature_departures(
                Cv, Cl, Hv, Hl, self._asplit_left, self._bsplit_left,
                N_stages, self.feed_enthalpies
            )
        return dTs
    
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
                    if j.specified_variable != 'Q':
                        specification_index.append(i)
                    missing.append(i)
                    continue
                bottom.phase = 'g'
                hv[i] = bottom.h
                bottom.phase = 'l'
            else:
                hv[i] = top.h
            if Li == 0:
                top.phase = 'l'
                hl[i] = top.h
                top.phase = 'g'
            else:
                hl[i] = bottom.h
            if j.specified_variable != 'Q':
                specification_index.append(i)
        if missing:
            neighbors = MESH.get_neighbors(missing=missing, size=N_stages)
            hv = MESH.fillmissing(neighbors, hv)
            hl = MESH.fillmissing(neighbors, hl)
        feed_enthalpies = self.feed_enthalpies
        if self.stage_reactions:
            feed_enthalpies = feed_enthalpies.copy()
            for i in self.stage_reactions:
                partition = partitions[i]
                feed_enthalpies[i] += partition.ins[0].Hf - sum([i.Hf for i in partition.outs])
        specification_index = np.array(specification_index, dtype=int)
        dB = MESH.phase_ratio_departures(
            L, V, hl, hv, 
            self._asplit_1, 
            self._asplit_left,
            self._bsplit_left,
            N_stages,
            specification_index,
            feed_enthalpies,
        )
        return dB
        
    def update_energy_balance_phase_ratio_departures(self):
        dBs = self.get_energy_balance_phase_ratio_departures()
        partitions = self.partitions
        # Bs = np.array([i.B for i in self.partitions])
        # Bs_new = dBs + Bs
        # index = np.argmin(Bs_new)
        # Bs_min = Bs_new[index]
        # if Bs_min < 0:
        #     f = - 0.99 * Bs[index] / dBs[index]
        #     print(f)
        #     breakpoint()
        #     dBs = f * dBs
        f = 1
        for i, dB in zip(partitions, dBs): 
            if dB >= 0: continue
            f = min(-0.99 * i.B / dB, f)
        dBs *= f
        for i, dB in zip(partitions, dBs):
            if i.specified_variable != 'B': 
                i.B += (1 - i.B_relaxation_factor) * dB
        if getattr(self, 'tracking', False):
            self._collect_variables('energy')
    
    def update_energy_balance_temperatures(self):
        dTs = self.get_energy_balance_temperature_departures()
        partitions = self.partitions
        for p, dT in zip(partitions, dTs):
            if p.specified_variable != 'T': 
                dT = (1 - p.T_relaxation_factor) * dT
                p.T += dT
                for i in p.outs: i.T += dT
        if getattr(self, 'tracking', False):
            self._collect_variables('energy')
        return dTs
    
    # %% Material balance convergence
    
    def _feed_flows_and_conversion(self):
        feed_flows = self.feed_flows.copy()
        partition = self.partitions
        index = self._update_index
        for i in self.stage_reactions: 
            p = partition[i]
            dmol = p.dmol
            for n, j in enumerate(index): feed_flows[i, n] += dmol[j]
        return feed_flows
    
    def set_flow_rates(self, bottom_flows, update_B=True):
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        feed_flows = self.feed_flows
        index = self._update_index
        if self.stage_reactions:
            feed_flows = self._feed_flows_and_conversion()
        f = PhasePartition.F_relaxation_factor
        if f and self.bottom_flows is not None:
            bottom_flows = f * self.bottom_flows + (1 - f) * bottom_flows
        self.bottom_flows = bottom_flows
        top_flows = MESH.top_flows_mass_balance(
            bottom_flows, feed_flows, self._asplit_left, self._bsplit_left, 
            self.N_stages
        )
        for i in range_stages:
            stage = stages[i]
            partition = stage.partition
            s_top, s_bot = partition.outs
            t = top_flows[i]
            mask = t < 0
            bulk_t = t.sum()
            if mask.any():
                t[mask] = 0
                t *= bulk_t / t.sum()
            b = bottom_flows[i]
            mask = b < 0
            bulk_b = b.sum()
            if mask.any():
                b[mask] = 0
                b *= bulk_b / b.sum()
            s_top.mol[index] = t
            s_bot.mol[index] = b
            for i in stage.splitters: i._run()
            if update_B and stage.specified_variable != 'B':
                stage.B = bulk_t / bulk_b
    
    def run_mass_balance(self):
        S = np.array([i.S for i in self.stages])
        feed_flows, *args = self._iter_args
        if self.stage_reactions:
            feed_flows = self._feed_flows_and_conversion()
        return MESH.bottom_flow_rates(S, feed_flows, *args)
       
    def update_mass_balance(self):
        self.set_flow_rates(self.run_mass_balance())
    
    def update_flow_rates(self, separation_factors, update_B=True):
        if separation_factors.min() < 0:
            S = separation_factors
            S_old = np.array([i.S for i in self._S_stages])
            # S * x + (1 - x) * S_old = 0
            # S * x + S_old - x * S_old = 0
            # S_old + x * (S - S_old) = 0
            # x = S_old / (S_old - S)
            x = S_old / (S_old - S)
            x = 0.1 * (x[(x < 1) & (x > 0)]).min()
            separation_factors = S * x + (1 - x) * S_old
        for stage, S in zip(self._S_stages, separation_factors): stage.S = S
        flows = self.run_mass_balance()
        self.set_flow_rates(flows, update_B)
        if getattr(self, 'tracking', False):
            self._collect_variables('material')
    
    # %% Phenomena convergence
    
    def update_bubble_point(self):
        stages = self.stages
        P = self.P
        if self.stage_reactions:
            self.update_liquid_holdup() # Finds liquid volume at each stage
            for stage in stages:
                partition = stage.partition
                partition._run_decoupled_KTvle(P=P)
                T = partition.T
                for i in (partition.outs + stage.outs): i.T = T
                if partition.reaction: 
                    partition._run_decoupled_reaction(P=P)
        else:
            for i in stages:
                partition = i.partition
                partition._run_decoupled_KTvle(P=P)
                T = partition.T
                for i in (partition.outs + i.outs): i.T = T
        if getattr(self, 'tracking', False):
            self._collect_variables('vle_phenomena')
    
    def update_pseudo_lle(self, separation_factors=None):
        stages = self.stages
        P = self.P
        if 'K' in self.partition_data:
            if separation_factors is not None: 
                self.update_flow_rates(separation_factors, update_B=False)
                for i in stages: i._update_separation_factors()
                if getattr(self, 'tracking', False):
                    self._collect_variables('material')
        else:
            stages = self._S_stages
            if separation_factors is None: separation_factors = np.array([i.S for i in stages])
            def psuedo_equilibrium(separation_factors):
                self.update_flow_rates(separation_factors, update_B=False)
                for n, i in enumerate(stages): 
                    i.partition._run_decoupled_Kgamma(P=P)
                    i._update_separation_factors()
                return np.array([i.S for i in stages])
            separation_factors = flx.fixed_point(
                psuedo_equilibrium, separation_factors, 
                xtol=self.S_tolerance,
                rtol=self.relative_S_tolerance,
                checkiter=False,
                checkconvergence=False,
            )
            self.update_flow_rates(separation_factors, update_B=True)
            for i in stages: i._update_separation_factors()
        for i in stages: 
            mixer = i.mixer
            partition = i.partition
            mixer.outs[0].mix_from(mixer.ins, energy_balance=False)
            partition._run_decoupled_B()
            i._update_separation_factors()
        if getattr(self, 'tracking', False):
            self._collect_variables('lle_phenomena')
    
    def update_pseudo_vle(self, separation_factors):
        P = self.P
        if 'K' in self.partition_data:
            if separation_factors is not None: 
                self.update_flow_rates(separation_factors)
                if getattr(self, 'tracking', False):
                    self._collect_variables('material')
        else:
            if separation_factors is None: separation_factors = np.array([i.S for i in self._S_stages])
            self.update_flow_rates(separation_factors, update_B=True)
            for i in self.stages: 
                i.partition._run_decoupled_Kfgas(P=P)
                i._update_separation_factors()
        if getattr(self, 'tracking', False):
            raise NotImplementedError('tracking sum rates decomposition not implemented')
        
