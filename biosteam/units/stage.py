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
from scipy.interpolate import LinearNDInterpolator, CubicSpline
from math import inf
from typing import Callable, Optional
from copy import copy
from scipy.optimize import fsolve, least_squares
# from cyipopt import minimize_ipopt
from scipy.optimize._numdiff import approx_derivative
from scipy.differentiate import jacobian
from scipy.spatial.distance import cdist
from collections import deque
from .. import Unit
from .design_tools import MESH
from numba import float64, int8, types, njit
from typing import NamedTuple, Iterable
import matplotlib.cm as cm
import matplotlib.colors as clr
import matplotlib.pyplot as plt
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

class IterationResult(NamedTuple):
    t: float #: Step size
    x: np.ndarray #: Point
    r: float #: Residual

# %% Inside-out tools

@jitdata
class SurrogateStage:
    _A: float
    _B: float
    _a: float64[:]
    _b: float64[:]
    _hV_ref: float
    _hL_ref: float
    _CV: float
    _CL: float
    _alpha: float64[:] 
    _beta: float64[:] #: K / (Kb * gamma) # For highly nonideal liquids
    _Kbmax: float
    _Kbmin: float
    
    def __init__( 
            self,
            x: float64[:], 
            x0: float64[:], 
            x1: float64[:],
            T: float, 
            T0: float, 
            T1: float,
            gamma: float64[:], 
            gamma0: float64[:], 
            gamma1: float64[:], 
            K: float64[:], 
            Kb: float64[:], 
            logK0: float64[:],
            logK1: float64[:],
            alpha: float64[:],
            beta: float64[:],
            w: float64[:],
            y: float64[:],
            hV: float, 
            hL: float, 
            CV: float, 
            CL: float, 
        ):
        b = (gamma1 - gamma0) / (x1 - x0)
        a = gamma - b * x
        Kb0 = np.exp((logK0 * w).sum())
        Kb1 = np.exp((logK1 * w).sum())
        B = np.log(Kb1 / Kb0) / (1/T1 - 1/T0)
        A = np.log(Kb) - B / T
        hV_ref = hV - CV * T
        hL_ref = hL - CL * T
        self._A = A
        self._B = B
        self._a = a
        self._b = b
        self._hV_ref = hV_ref
        self._hL_ref = hL_ref
        self._CV = CV
        self._CL = CL
        self._alpha = alpha
        self._beta = beta
        if Kb0 < Kb1: 
            self._Kbmin = Kb0
            self._Kbmax = Kb1
        else:
            self._Kbmin = Kb1
            self._Kbmax = Kb0
    
    def gamma(self, x: float64[:], T: float) -> float64[:]:
        return self._a + self._b * x 
    
    def Kb(self, T: float) -> float:
        return np.exp(self._A + self._B / T)

    def T(self, Kb: float) -> float:
        return self._B / (np.log(Kb) - self._A)
    
    def bubble_point(self, x: float64[:]) -> types.UniTuple(float, 2):
        Kbmin = self._Kbmin
        Kbmax = self._Kbmax
        beta = self._beta
        fgamma = self.gamma
        fT = self.T
        zero = x < 0
        if zero.any():
            x[zero] = 0
            x /= x.sum()
        
        # Initial guess without gamma
        Kb = 1 / (self._alpha * x).sum()
        if Kb < Kbmin: Kb = Kbmin
        if Kb > Kbmax: Kb = Kbmax
        T0 = fT(Kb)
        
        # Guess with gamma
        alpha = fgamma(x, T0) * beta
        Kb = 1 / (alpha * x).sum()
        if Kb < Kbmin: Kb = Kbmin
        if Kb > Kbmax: Kb = Kbmax
        T1 = g0 = self.T(Kb)
        
        # Solve by accelerated Wegstein iteration
        for _ in range(100):
            dT = T1 - T0
            alpha = fgamma(x, T0) * beta
            Kb = 1 / (alpha * x).sum()
            if Kb < Kbmin: Kb = Kbmin
            if Kb > Kbmax: Kb = Kbmax
            g1 = fT(Kb)
            error = np.abs(g1 - T1)
            if error < 1e-9: 
                T = g1
                break
            T0 = T1
            denominator = dT - g1 + g0
            if abs(denominator) > 1e-16:
                w = (dT / denominator)
                if not w < 0: w = w ** 0.5
                T1 = w*g1 + (1.-w)*T1
            else:
                T1 = g1
            g0 = g1
        else:
            T = g1
        return Kb, T

    def hV(self, T: float) -> float:
        return self._hV_ref + self._CV * T
    
    def hL(self, T: float) -> float:
        return self._hL_ref + self._CL * T

@jitdata
class SurrogateColumn:
    stages: types.List(SurrogateStage.class_type.instance_type, reflected=True)
    N_stages: int
    N_chemicals: int
    specified_variables: str
    specified_values: float64[:]
    neg_asplit: float64[:]
    neg_bsplit: float64[:]
    top_split: float64[:]
    bottom_split: float64[:]
    feed_and_invariable_enthalpies: float64[:] 
    feed_flows: float64[:, :]
    total_feed_flows: float64[:]
    bulk_feed: float
    alpha: float64[:, :]
    Kb: float64[:]
    point_shape: types.UniTuple(int8, 2)
    
    def __init__(
            self,
            N_stages: int, 
            N_chemicals: int,
            T: float, 
            x: float64[:, :], 
            y: float64[:, :], 
            gamma: float64[:, :], 
            K: float64[:, :], 
            dlogK_dTinv: float64[:, :],
            hV: float64[:],
            hL: float64[:], 
            CV: float64[:],
            CL: float64[:],
            specified_variables: str,
            specified_values: float64[:],
            neg_asplit: float64[:],
            neg_bsplit: float64[:],
            top_split: float64[:],
            bottom_split: float64[:],
            feed_and_invariable_enthalpies: float64[:],
            feed_flows: float64[:, :],
            total_feed_flows: float64[:],
            bulk_feed: float,
            point_shape: types.UniTuple(int8, 2),
        ):
        logK = np.log(K + 1e-34)
        w = y * dlogK_dTinv
        w /= np.expand_dims(w.sum(axis=1), -1)
        Kb = np.exp((logK * w).sum(axis=1))
        alpha = K / np.expand_dims(Kb, -1)
        beta = alpha / gamma
        last = N_stages - 1
        stages = []
        for i in range(N_stages):
            if i == 0:
                i0 = i
                i1 = i + 1
            elif i == last:
                i0 = i - 1
                i1 = i
            else:
                i0 = i - 1
                i1 = i + 1
            stages.append(
                SurrogateStage(
                    x[i], x[i0], x[i1],
                    T[i], T[i0], T[i1],
                    gamma[i], gamma[i0], gamma[i1],
                    K[i], Kb[i], logK[i0], logK[i1], 
                    alpha[i], beta[i], w[i], y[i], 
                    hV[i], hL[i], CV[i], CL[i],
                )
            )
        self.stages = stages 
        self.N_stages = N_stages 
        self.N_chemicals = N_chemicals
        self.specified_variables = specified_variables
        self.specified_values = specified_values
        self.neg_asplit = neg_asplit
        self.neg_bsplit = neg_bsplit
        self.top_split = top_split
        self.bottom_split = bottom_split
        self.feed_and_invariable_enthalpies = feed_and_invariable_enthalpies
        self.feed_flows = feed_flows
        self.total_feed_flows = total_feed_flows
        self.bulk_feed = bulk_feed
        self.alpha = alpha
        self.Kb = Kb
        self.point_shape = point_shape
    
    def Sb_to_point(self, Sb: float64[:]) -> float64[:, :]:
        stages = self.stages
        N_stages = self.N_stages
        S = self.alpha * np.expand_dims(Sb, -1)
        xL = MESH.bottom_flow_rates(
            S, 
            self.feed_flows, 
            self.neg_asplit, 
            self.neg_bsplit,
            N_stages,
        )
        L = xL.sum(axis=1)
        x = xL / np.expand_dims(L, -1)
        yV = xL * S
        T = np.zeros(N_stages)
        for i, stage in enumerate(stages):
            _, Ti = stage.bubble_point(x[i])
            T[i] = Ti
        n = self.N_chemicals
        point = np.zeros(self.point_shape)
        point[:, :n] = yV
        point[:, n] = T
        point[:, -n:] = xL
        return point
    
    def residuals(self, logSb1: float64[:]) -> float64[:]:
        stages = self.stages
        N_stages = self.N_stages
        Sb = np.exp(logSb1) - 1
        S = self.alpha * np.expand_dims(Sb, -1)
        xL = MESH.bottom_flow_rates(
            S, 
            self.feed_flows, 
            self.neg_asplit, 
            self.neg_bsplit,
            N_stages,
        )
        L = xL.sum(axis=1)
        x = xL / np.expand_dims(L, -1)
        yV = xL * S
        V = yV.sum(axis=1)
        T = np.zeros(N_stages)
        hV = T.copy()
        hL = T.copy()
        residuals = T.copy()
        for i, stage in enumerate(stages):
            Kb, Ti = stage.bubble_point(x[i])
            T[i] = Ti
            hV[i] = stage.hV(Ti)
            hL[i] = stage.hL(Ti)
        HV = hV * V
        HL = hL * L
        variables = self.specified_variables
        values = self.specified_values
        H_feed = self.feed_and_invariable_enthalpies
        top_split = self.top_split
        bottom_split = self.bottom_split
        for i, stage in enumerate(stages):
            var = variables[i]
            value = values[i]
            if var == 'Q':
                H_out = HV[i] + HL[i]
                H_in = H_feed[i]
                i0 = i - 1
                if i0 != -1: 
                    H_in += (1 - bottom_split[i0]) * HL[i0]
                i1 = i + 1
                if i1 != N_stages: 
                    H_in += (1 - top_split[i1]) * HV[i1]
                residuals[i] = H_out - H_in - value
            elif var == 'T':
                residuals[i] = value - T[i]
            elif var == 'B':
                residuals[i] = V[i] - L[i] * value
            elif var == 'F':
                residuals[i] = value * self.bulk_feed - L[i]
            else:
                raise RuntimeError('unknown specification')
        return residuals
    
    def phenomena_iter(self, Sb: float64[:]) -> float64[:]:
        N_stages = self.N_stages
        S = self.alpha * np.expand_dims(Sb, -1)
        xL = MESH.bottom_flow_rates(
            S, 
            self.feed_flows, 
            self.neg_asplit, 
            self.neg_bsplit,
            N_stages,
        )
        L = xL.sum(axis=1)
        x = xL / np.expand_dims(L, -1)
        stages = self.stages
        hL = np.zeros(N_stages)
        hV = hL.copy() 
        Kb = hL.copy()
        T = hL.copy()
        for i in range(N_stages):
            stage = stages[i]
            Kbi, Ti = stage.bubble_point(x[i])
            Kb[i] = Kbi
            T[i] = Ti
            hL[i] = stage.hL(Ti)
            hV[i] = stage.hV(Ti)
        V, L = MESH.bulk_vapor_and_liquid_flow_rates(
            hL, hV,
            self.neg_asplit, self.neg_bsplit, 
            self.top_split, self.bottom_split, 
            N_stages, self.feed_and_invariable_enthalpies, 
            self.total_feed_flows,
            self.specified_variables,
            self.specified_values,
            self.bulk_feed,
        )
        return Kb * V / L

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


@njit(types.UniTuple(float64[:, :, :], 3)(types.List(JacobianData.class_type.instance_type, reflected=True), int8, int8, int8))
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
    
    def _init(self, T=None, P=None, Q=None, phase=None):
        self.specify_variables(T, P, Q)
        self.T = T
        self.Q = Q
        self.P = P
        self.phase = phase
        
    def specify_variables(self, T=None, P=None, Q=None):
        if T is not None:
            self.specified_variable = 'T'
        elif Q is not None:
            self.specified_variable = 'Q'
        else:
            self.specified_variable = 'Q'
            Q = 0
        self.T = T
        self.P = P
        self.Q = Q
        
    def _run(self):
        outlet = self.outs[0]
        outlet.mix_from(self.ins, energy_balance=False)
        if self.P is not None: outlet.P = self.P
        if self.phase is None: 
            outlet.phase = self.ins[0].phase
        else:
            outlet.phase = self.phase
        if self.specified_variable == 'Q':
            outlet.H = sum([i.H for i in self.ins], self.Q)
        else:
            outlet.T = self.T

    def _update_energy_coefficient(self, stream, coefficients):
        if self.specified_variable == 'Q':
            C = stream.C
            coefficients[stream, 'T'] = -C
            return C * stream.T
        return 0
    
    def _create_bulk_balance_equations(self):
        fresh_inlets, process_inlets, equations = self._begin_bulk_equations()
        outlet = self.outs[0]
        coeff = {(outlet, 'F_mol'): 1}
        for i in process_inlets: coeff[i, 'F_mol'] = -1
        return [(coeff, sum([i.F_mol for i in fresh_inlets]))]
    
    def _create_energy_balance_equations(self):
        if self.specified_variable == 'Q':
            fresh_inlets, process_inlets, equations = self._begin_energy_equations()
            outlet = self.outs[0]
            C = outlet.C
            coeff = {(self, 'T'): C,
                     (outlet, 'F_mol'): outlet.h}
            Q = self.Q + outlet.T * C + sum([i.H for i in fresh_inlets])
            for i in process_inlets: 
                coeff[i, 'F_mol'] = -i.h
                Q -= i._update_energy_coefficient(coeff)
            return [(coeff, Q)]
        else:
            return []
        
    def _create_material_balance_equations(self, composition_sensitive):
        outlet = self.outs[0]
        fresh_inlets, process_inlets, equations = self._begin_material_equations(composition_sensitive)
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

    def _update_variable(self, variable, value):
        self.outs[0].T = value
        
    def _update_nonlinearities(self): pass
    
    @property
    def equation_node_names(self): 
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
        if self.specified_variable == 'Q':
            self.energy_balance_node.set_equations(
                inputs=(
                    self.T_node, 
                    *[i.T_node for i in (*self.ins, *self.outs)],
                    *[i.F_node for i in (*self.ins, *self.outs)],
                    *[j for i in self.ins if (j:=i.E_node)]
                ),
                outputs=[j for i in self.outs if (j:=i.E_node)],
            )
        else:
            self.energy_balance_node.set_equations(
                inputs=[i.F_node for i in (*self.ins, *self.outs)],
            )
    
    @property
    def T_node(self):
        if hasattr(self, '_T_node'): return self._T_node
        self._T_node = var = VariableNode(f"{self.node_tag}.T", lambda: self.T)
        return var 
        
    @property
    def E_node(self):
        if self.specified_variable == 'Q':
            return None
        else:
            return self.T_node
    

class ReactivePhaseStage(bst.Unit): # Does not include VLE
    _N_outs = _N_ins = 1
    _ins_size_is_fixed = False
    
    @property
    def equation_node_names(self): 
        return (
            'overall_material_balance_node', 
            'reaction_phenomenode',
            'energy_balance_node',
        )
    
    def _init(self, reaction, T=None, P=None, Q=None, phase=None):
        self.specify_variables(T, P, Q)
        self.reaction = reaction
        self.phase = phase
        
    def specify_variables(self, T=None, P=None, Q=None):
        if T is not None:
            self.specified_variable = 'T'
        elif Q is not None:
            self.specified_variable = 'Q'
        else:
            self.specified_variable = 'Q'
            Q = 0
        self.T = T
        self.P = P
        self.Q = Q
        
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
        fresh_inlets, process_inlets, equations = self._begin_material_equations(composition_sensitive)
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
    
    def _update_variable(self, variable, value):
        self.outs[0].T = value
    
    def _update_energy_coefficient(self, stream, coefficients):
        if self.specified_variable == 'Q':
            C = stream.C
            coefficients[stream, 'T'] = -C
            return C * stream.T
        return 0
    
    def _create_bulk_balance_equations(self):
        fresh_inlets, process_inlets, equations = self._begin_bulk_equations()
        outlet = self.outs[0]
        coeff = {(outlet, 'F_mol'): 1}
        for i in process_inlets: coeff[i, 'F_mol'] = -1
        return [(coeff, self.dmol.sum() + sum([i.F_mol for i in fresh_inlets]))]
    
    def _create_energy_balance_equations(self):
        if self.specified_variable == 'Q':
            fresh_inlets, process_inlets, equations = self._begin_energy_equations()
            outlet = self.outs[0]
            C = outlet.C
            coeff = {(self, 'T'): C,
                     (outlet, 'F_mol'): outlet.h + outlet.hf}
            Q = self.Q + C * self.T + sum([i.H + i.Hf for i in fresh_inlets])
            for i in process_inlets: 
                coeff[i, 'F_mol'] = -i.h - i.hf
                Q -= i._update_energy_coefficient(coeff)
            return [(coeff, Q)]
        else:
            return []
    
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
        if self.specified_variable == 'Q':
            self.energy_balance_node.set_equations(
                inputs=(
                    self.T_node, 
                    *[i.T_node for i in (*self.ins, *self.outs)],
                    *[i.F_node for i in (*self.ins, *self.outs)],
                    *[j for i in self.ins if (j:=i.E_node)]
                ),
                outputs=[j for i in self.outs if (j:=i.E_node)],
            )
        else:
            self.energy_balance_node.set_equations(
                inputs=[i.F_node for i in (*self.ins, *self.outs)],
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
        if self.specified_variable == 'Q':
            return self.E_node
        else:
            return None


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
    
    def _update_auxiliaries(self):
        for i in self.splitters: i.ins[0].mix_from(i.outs, energy_balance=False)
        self.mixer.outs[0].mix_from(self.ins, energy_balance=False)
    
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
    def top(self):
        return self.outs[0]
    @property
    def bottom(self):
        return self.outs[1]
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
    def F(self):
        return self.partition.F
    @F.setter
    def F(self, F):
        self.partition.F = F
    
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
    def x(self):
        return self.partition.x
    @x.setter
    def x(self, x):
        self.partition.x = x
        
    @property
    def y(self):
        return self.partition.y
    @y.setter
    def y(self, y):
        self.partition.y = y
    
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
        
    def _update_separation_factors(self, f=None):
        if self.B == inf:
            self.S = np.ones(len(self.partition.IDs)) * np.inf
        elif self.B == 0: 
            self.S = np.zeros(len(self.partition.IDs))
        else:
            K = self.K
            S = K * self.B
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
    def K_model(self):
        try:
            K_model = self._K_model
        except:
            if 'K' in self.partition.partition_data:
                self._K_model = K_model = lambda y, x, T, P: self.partition.partition_data['K']
            else:
                self._K_model = K_model = tmo.equilibrium.PartitionCoefficients(''.join(self.phases), self.chemicals[self.partition.IDs], self.thermo)
        return K_model
    
    @property
    def _equilibrium_residuals_vectorized(self):
        try:
            return self._vectorized_equilibrium_residuals
        except:
            K_model = self.K_model
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
                    error = K * mol_bottom * bulk_top / bulk_bottom - mol_top
                elif self.B == 0 or bulk_top:
                    error = -0.1 * mol_top
                else:
                    error = 0.1 * mol_bottom
                # if 'g' in self.phases:
                #     try:
                #         Tmax = self._Tmax
                #         Tmin = self._Tmin
                #     except:
                #         chemicals = self.chemicals[self.partition.IDs]
                #         Ts = [i.Tsat(self.P) for i in chemicals]
                #         Tmax = self._Tmax = 1.2 * max(Ts)
                #         Tmin = self._Tmin = 0.8 * min(Ts)
                #     if T > Tmax:
                #         error[error > 0] += 1e3 * (T - Tmax)
                #         error[error < 0] -= 1e3 * (T - Tmax)
                #     elif T < Tmin:
                #         error[error > 0] += 1e3 * (Tmin - T)
                #         error[error < 0] -= 1e3 * (Tmin - T)
                return error
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
            K_model = self.K_model
            K = K_model(y, x, T, self.P)
            error = K * mol_bottom * bulk_top / bulk_bottom - mol_top
        elif self.B == 0 or bulk_top:
            error = -0.1 * mol_top
        else:
            error = 0.1 * mol_bottom
        return error
    
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
            F = center.bottom.mol.sum()
            return self.F * self._bulk_feed - F
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
        index = self._eq_index
        N_chemicals = self._N_chemicals
        top, bottom = self.partition.outs
        top.mol[index] = mol_top = x[:N_chemicals]
        if self.specified_variable != 'T':
            self.T = x[N_chemicals]
        bottom.mol[index] = mol_bot = x[-N_chemicals:]
        if self.specified_variable != 'B':
            bulk_bot = mol_bot.sum()
            if bulk_bot:
                self.B = mol_top.sum() / bulk_bot
            else:
                self.B = float('inf')
        mol_bot[mol_bot == 0] = 1-16
        self.S = mol_top / mol_bot
        if self.B == 0: 
            self.K = np.zeros(N_chemicals)
        else:
            self.K = self.S / self.B
        for i in self.splitters: i._run()

    # %% Phenomena-based simulation
    
    @property
    def composition_sensitive(self):
        return self.phases == ('L', 'l')
    
    def _update_energy_coefficient(self, stream, coefficients):
        if self.specified_variable == 'Q' and self.phases == ('L', 'l'):
            C = stream.C
            coefficients[self, 'T'] = -C
            H_ref = C * self.T
        else:
            H_ref = 0
        return H_ref
    
    def _create_energy_balance_equations(self):
        if self.specified_variable != 'Q': return []
        fresh_inlets, process_inlets, equations = self._begin_energy_equations()
        Q = self.Q + sum([i.H for i in fresh_inlets])
        if self.reaction:
            coeff = {(i, 'F_mol'): i.h + i.hf for i in self.outs}
        else:
            coeff = {(i, 'F_mol'): i.h for i in self.outs}
        if self.phases == ('L', 'l'):
            coeff[self, 'T'] = C = sum([i.C for i in self.outs])
            Q += C * self.T
        for i in process_inlets: 
            if i.isempty(): continue
            if self.reaction:
                coeff[i, 'F_mol'] = -i.h - i.hf
            else:
                coeff[i, 'F_mol'] = -i.h
            Q -= i._update_energy_coefficient(coeff)
        return [(coeff, Q)]
    
    def _create_bulk_balance_equations(self):
        top_split = self.top_split
        bottom_split = self.bottom_split
        fresh_inlets, process_inlets, equations = self._begin_bulk_equations()
        top, bottom, *_ = self.outs
        top_side_draw = self.top_side_draw
        bottom_side_draw = self.bottom_side_draw

        # # Overall flows
        eq_overall = {}
        reaction = self.reaction
        if reaction: # Reactive liquid
            predetermined_flow = sum([i.F_mol for i in fresh_inlets])
            rhs = predetermined_flow + self.partition.dmol.sum()
            for i in self.outs: eq_overall[i, 'F_mol'] = 1
            for i in process_inlets: eq_overall[i, 'F_mol'] = -1
            equations.append(
                (eq_overall, rhs)
            )
        else:
            for i in self.outs: eq_overall[i, 'F_mol'] = 1
            for i in process_inlets:
                if i in eq_overall: del eq_overall[i, 'F_mol']
                else: eq_overall[i, 'F_mol'] = -1
            equations.append(
                (eq_overall, sum([i.F_mol for i in fresh_inlets]))
            )
        if self.specified_variable != 'Q' or self.phases == ('L', 'l'):
            eq_outs = {}
            # self._update_auxiliaries()
            # partition = self.partition
            # top, bottom = self.partition.outs
            # IDs = partition.IDs
            # B = top.imol[IDs].sum() / bottom.imol[IDs].sum()
            B = self.B
            if top_split == 1:
                if bottom_split == 1:
                    eq_outs[top_side_draw, 'F_mol'] = -1
                    eq_outs[bottom_side_draw, 'F_mol'] = B
                else:
                    eq_outs[top_side_draw, 'F_mol'] = bottom_split - 1
                    eq_outs[bottom, 'F_mol'] = B
            elif bottom_split == 1:
                eq_outs[top, 'F_mol'] = -1
                eq_outs[bottom_side_draw, 'F_mol'] = B * (1 - top_split) 
            else:
                eq_outs[top, 'F_mol'] = bottom_split - 1 
                eq_outs[bottom, 'F_mol'] = B * (1 - top_split) 
            equations.append(
                (eq_outs, 0)
            )
            # Top split flows
            if top_side_draw:
                if top_split == 1:
                    eq_top_split = {
                        (top, 'F_mol'): 1,
                    }
                else:
                    eq_top_split = {
                        (top_side_draw, 'F_mol'): 1,
                        (top, 'F_mol'): -top_split / (1 - top_split),
                    }
                equations.append(
                    (eq_top_split, 0)
                )
            # Bottom split flows
            if bottom_side_draw:
                if bottom_split == 1:
                    eq_bottom_split = {
                        (bottom, 'F_mol'): 1,
                    }
                else:
                    eq_bottom_split = {
                        (bottom_side_draw, 'F_mol'): 1,
                        (bottom, 'F_mol'): -bottom_split / (1 - bottom_split),
                    }
                equations.append(
                    (eq_bottom_split, 0)
                )
        return equations
    
    def _create_material_balance_equations(self, composition_sensitive):
        self._update_separation_factors()
        partition = self.partition
        chemicals = self.chemicals
        pIDs = partition.IDs
        IDs = chemicals.IDs
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
        fresh_inlets, process_inlets, equations = self._begin_material_equations(composition_sensitive)
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
    
    def _reset_bulk_variable(self):
        self._update_auxiliaries()
        if self.specified_variable != 'B' and self.phases != ('L', 'l'):
            partition = self.partition
            top, bottom = self.partition.outs
            IDs = partition.IDs
            F_mol_top = top.imol[IDs].sum()
            F_mol_bottom = bottom.imol[IDs].sum()
            if F_mol_bottom:
                partition.B = F_mol_top / F_mol_bottom
            elif F_mol_top:
                partition.B = 1e3
            else:
                pass
            
    def _update_variable(self, variable, value):
        if variable == 'T':
            self.T = value
            for i in self.outs: i.T = value
        else:
            raise ValueError('cannot update variable')
            
    def _update_composition_parameters(self):
        partition = self.partition
        data = partition.partition_data
        if data and 'K' in data: return
        partition._run_decoupled_Kgamma()
    
    def _update_net_flow_parameters(self):
        mixer = self.mixer
        mixer.outs[0].mix_from(mixer.ins, energy_balance=False)
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
        balances = (
            'overall_material_balance_node', 
            'separation_material_balance_node',
            'energy_balance_node',
        )
        if self.phases == ('g', 'l'):
            phenomenode = 'vle_phenomenode'
        else: # Assume LLE
            phenomenode = 'lle_phenomenode'
        return (
            *balances,
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
        elif  (partition_data and 'K' in partition_data):
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
        if self.phases == ('g', 'l'):
            return self.Phi_node
        else:
            return self.T_node
    
    
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
        self.dmol = SparseVector.from_size(self.chemicals.size)
        for i, j in zip(self.outs, self.phases): i.phase = j 
        
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
                    fname = name + '_relaxation_factor'
                    if hasattr(self, fname):
                        f = getattr(self, fname)
                        g = 1. - f 
                        for i, j in enumerate(index):
                            last[j] = g * array[i] + f * last[j]
                    else:
                        for i, j in enumerate(index):
                            last[j] = array[i]
            else:
                self.IDs = IDs
                index = [IDs.index(i) for i in IDs_last]
                for name, array in kwargs.items():
                    last = getattr(self, name)
                    new = array.copy()
                    setattr(self, name, new)
                    for i, j in enumerate(index): new[i] = last[j]
                    fname = name + '_relaxation_factor'
                    if hasattr(self, fname):
                        f = getattr(self, fname)
                        g = 1. - f 
                        for i, j in enumerate(index):
                            new[j] = g * array[i] + f * new[j]
                    else:
                        for i, j in enumerate(index):
                            new[j] = array[i]
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
            self._set_arrays(IDs, K=K_new, x=p.x, y=p.y)
    
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
        var = self.specified_variable
        if var == 'F':
            bottom = self.F * self._bulk_feed
            feed_mol = ms.mol.sum()
            self.B = feed_mol  / bottom - 1
            var = 'B'
        kwargs = {
            var: getattr(self, var),
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
            if not L_total: x_mol = y_mol
            K_new = y_mol / x_mol
        else:
            y_mol = x_mol
            K_new = np.ones(len(index)) * 1e-16
        if 'B' not in kwargs:
            if not L_total:
                self.B = inf
            else:
                self.B = V_total / L_total
        self.T = ms.T
        self._set_arrays(IDs, K=K_new, x=x_mol, y=y_mol)
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
    _line_search = False # Experimental feature.
    _tracked_points = None # For numerical/convergence analysis.
    minimum_residual_reduction = 0.25 # Minimum fractional reduction in residual for simulation.
    iteration_memory = 10 # Length of recorded iterations.
    inside_maxiter = 100
    default_max_attempts = 5
    default_maxiter = 100
    default_optimize_result = True
    default_tolerance = 1e-5
    default_relative_tolerance = 1e-4
    default_algorithms = ('inside out', 'phenomena', 'simultaneous correction', 'sequential modular')
    default_inside_loop_algorithm = 'simultaneous correction'
    decomposition_algorithms = {
        'phenomena', 'inside out', 'phenomena modular', 'sequential modular',
    }
    available_algorithms = {
        *decomposition_algorithms, 
        'simultaneous correction',
    }
    default_methods = {
        'phenomena': 'fixed-point',
        'phenomena modular': 'fixed-point',
        'sequential modular': 'fixed-point',
        'inside out': 'fixed-point',
        'simultaneous correction': 'hybr', # Alternatively 'trf'
    }
    method_options = {
        'fixed-point': {},
        'wegstein': {'lb': 1, 'ub': 4, 'exp': 0.5}
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
            algorithms=None,
            methods=None,
            maxiter=None,
            max_attempts=None,
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
        if 'K' in partition_data:
            K = partition_data['K']
            partition_stage = K.ndim == 2
        else:
            partition_stage = False
        self.N_stages = N_stages
        if not isinstance(P, Iterable): P = [P] * N_stages 
        self.multi_stream = tmo.MultiStream(None, P=P[N_stages // 2], phases=phases, thermo=self.thermo)
        self.P = np.array(P)
        self.T = T
        self.phases = phases = self.multi_stream.phases # Corrected order
        self._has_vle = 'g' in phases
        self._has_lle = 'L' in phases
        self._top_split = top_splits = np.zeros(N_stages)
        self._bottom_split = bottom_splits = np.zeros(N_stages)
        self._convergence_analysis_mode = False
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
                    outs.append(next(bsd_iter))
                    bottom_split = bottom_side_draws[i]
                    bottom_splits[i] = bottom_split
                else: 
                    bottom_split = 0
                if partition_stage:
                    pd = partition_data.copy()
                    pd['K'] = K[i]
                else:
                    pd = partition_data
                new_stage = self.auxiliary(
                    'stages', StageEquilibrium, phases=phases,
                    ins=feed,
                    outs=outs,
                    partition_data=pd,
                    top_split=top_split,
                    bottom_split=bottom_split,
                    P=P[i]
                )
                if last_stage:
                    last_stage.add_feed(new_stage-0)
                last_stage = new_stage
            for feed, stage in zip(self.ins, feed_stages):
                stages[stage].add_feed(self.auxlet(feed))  
            #: dict[int, tuple(str, float)] Specifications for VLE by stage
            self.stage_specifications = stage_specifications
            for i, (name, value) in stage_specifications.items():
                stages[i].specify_variables(**{name: value}, P=P[i])
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
        self._asplit = 1 - top_splits
        self._bsplit = 1 - bottom_splits
        self._neg_asplit = top_splits - 1
        self._neg_bsplit = bottom_splits - 1
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
        
        #: [bool] Optimize final result.
        self.optimize_result = self.default_optimize_result

        #: [float] Absolute molar flow and temperature tolerance
        self.tolerance = self.default_tolerance

        #: [float] Relative molar flow and temperature tolerance
        self.relative_tolerance = self.default_relative_tolerance
        
        self.inside_loop_algorithm = self.default_inside_loop_algorithm
        
        self.use_cache = True if use_cache else False
        
        self.collapsed_init = collapsed_init
        
        if algorithms is None:
            self.algorithms = self.default_algorithms
        elif isinstance(algorithms, str):
            self.algorithms = (algorithms,)
        else:
            self.algorithms = algorithms        
            
        if methods is None:
            self.methods = [self.default_methods[i] for i in self.algorithms]
        elif isinstance(methods, str):
            self.methods = len(self.algorithms) * (methods,)
        else:
            self.methods = methods
        
        self.vle_decomposition = vle_decomposition
    
    
    # %% Optimization
    
    def _get_point(self, x=None):
        if x is None: x = np.zeros(self._point_shape)
        for i, stage in enumerate(self.stages): 
            stage._get_point(x[i])
        return x
    
    def _set_point(self, x):
        for i, stage in enumerate(self.stages): 
            stage._set_point(x[i])
    
    def _objective(self, x):
        return self._net_residual(self._residuals(x))
    
    def _net_residual(self, residuals):
        return (residuals * residuals).sum()
    
    def _residuals(self, x):
        x[x < 0] = 0
        N_chemicals = self._N_chemicals
        N_stages, N_variables = x.shape
        stages = self.stages
        residuals = np.zeros(N_variables * N_stages) # H, Mi, Ei
        H_magnitude = self._H_magnitude
        stage_data = [
            stage._stage_data(xi, H_feed, mol_feed, H_magnitude)
            for xi, stage, H_feed, mol_feed
            in zip(x, stages, self._feed_and_invariable_enthalpies, self.feed_flows)
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
        i = 1
        stage = stages[i]
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
        for i, (xi, stage, H_feed, mol_feed) in enumerate(zip(x, stages, self._feed_and_invariable_enthalpies, self.feed_flows)):
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
    
    def _run_simultaneous_correction(self, x): 
        # Newton-Raphson with an exact line search.
        # fsolve (which uses trust-region) is better.
        # This method marely serves as a baseline to perform numerical analysis.
        residuals = self._residuals(x)
        A, B, C = self._jacobian(x)
        correction = -MESH.solve_block_tridiagonal_matrix(A, B, C, residuals)
        
        # N_chemicals = self._N_chemicals
        # N_stages, N_variables = x.shape
        # jacobian_data = []
        # stage_data = []       
        # stages = self.stages
        # H_magnitude = self._H_magnitude
        # residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
        # for i, (xi, stage, H_feed, mol_feed) in enumerate(zip(x, stages, self.feed_enthalpies, self.feed_flows)):
        #     JD, SD = stage._simultaneous_correction_data(xi, H_feed, mol_feed, H_magnitude)
        #     jacobian_data.append(JD)
        #     stage_data.append(SD)
        # center = stage_data[0]
        # lower = stage_data[1]
        # stage = stages[0]
        # H_index = 0
        # M_slice = slice(1, N_chemicals + 1)
        # E_slice = slice(N_chemicals + 1, None)
        # residuals = np.zeros([N_stages, N_variables]) # H, Mi, Ei
        # residuals[0, H_index] = stage._energy_balance_residual(None, center, lower)
        # residuals[0, M_slice] = stage._material_balance_residuals(None, center, lower)
        # stage = stages[1]
        # for i in range(2, N_stages): 
        #     upper = center
        #     center = lower
        #     lower = stage_data[i]
        #     ilast = i-1
        #     residuals[ilast, H_index] = stage._energy_balance_residual(upper, center, lower)
        #     residuals[ilast, M_slice] = stage._material_balance_residuals(upper, center, lower)
        #     residuals[ilast, E_slice] = stage._equilibrium_residuals(center)
        #     stage = stages[i]
        # upper = center
        # center = lower
        # residuals[i, H_index] = stage._energy_balance_residual(upper, center, None)
        # residuals[i, M_slice] = stage._material_balance_residuals(upper, center, None)
        # residuals[i, E_slice] = stage._equilibrium_residuals(center)
        # A, B, C = jacobian_blocks(jacobian_data, N_stages, N_chemicals, N_variables)
        # correction = -MESH.solve_block_tridiagonal_matrix(A, B, C, residuals)
        try:
            result = flx.inexact_line_search(
                self._objective, x, correction, t0=1e-3, t1=0.1,
            )
        except:
            return x + 1e-3 * correction
        else:
            return result.x
    
    # %% Decoupled phenomena equation oriented simulation
    
    def _update_auxiliaries(self):
        for i in self.stages: i._update_auxiliaries()
    
    @property
    def composition_sensitive(self):
        return self._has_lle
    
    def _update_composition_parameters(self):
        for i in self.partitions: 
            if 'K' in i.partition_data: continue
            i._run_decoupled_Kgamma()
    
    def _update_net_flow_parameters(self):
        for i in self.stages: i._update_net_flow_parameters()
    
    def _update_nonlinearities(self):
        if self._has_vle:
            for i in self.stages: i._update_nonlinearities()
        elif self._has_lle:
            pass
            # self.update_pseudo_lle()
    
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
        x = self.hot_start()
        algorithms = self.algorithms
        methods = self.methods
        optimize_result = self.optimize_result
        f = self._iter
        analysis_mode = self._convergence_analysis_mode
        maxiter = self.maxiter
        xtol = 100 * self.tolerance if optimize_result else self.tolerance
        rtol = 100 * self.relative_tolerance if optimize_result else self.relative_tolerance
        self.iter = 0
        for n in range(self.max_attempts):
            self.attempt = n
            for algorithm, method in zip(algorithms, methods):
                if algorithm == 'simultaneous correction':
                    x = self._simultaneous_correction(x, method)
                    if self._best_result.r < 1e-3:
                        optimize_result = False
                        break
                    continue
                if method == 'fixed-point':
                    solver = flx.fixed_point
                elif method == 'wegstein':
                    solver = flx.wegstein
                else:
                    raise ValueError(f'invalid method {method!r}')
                if self._has_lle and algorithm == 'inside out': continue
                if analysis_mode:
                    self._tracked_algorithms.append(
                        (self.iter + 1, algorithm)
                    )
                try: x = solver(f, x, maxiter=maxiter, xtol=xtol, rtol=rtol, args=(algorithm,))
                except:
                    self._mean_residual = np.inf
                    if self._best_result.x is not None: x = self._best_result.x
                    maxiter = self.maxiter - self.iter
                    if maxiter <= 0: break
                else: 
                    x = self._best_result.x 
                    break
            else: continue
            break
        if optimize_result: x = self._simultaneous_correction(x, 'hybr')
        self._set_point(x)
        # Last simulation to force mass balance
        self.update_mass_balance()
    
    def _new_point(self, x1=None):
        record = self._iteration_record
        if x1 is None: x1 = self._get_point()
        if self._line_search:
            t0, x0, r0 = record[0]
            correction = x1 - x0
            result = flx.inexact_line_search(
                self._objective, x0, correction, 
                fx=r0, t0=0.8, t1=1.2, tguess=t0
            )
        else:
            result = IterationResult(1, x1, self._objective(x1))
        x1 = result.x
        x1[x1 < 0] = 0
        record.rotate()
        record[0] = result
        residuals = np.array([i.r for i in record])
        mean = np.mean(residuals)
        mrr = self.minimum_residual_reduction
        if self._best_result.r > result.r: self._best_result = result
        if mean > self._mean_residual * (1 - mrr):
            raise RuntimeError('residual error is not decreasing sufficiently')
        else:
            self._mean_residual = mean
            return x1
    
    def _iter(self, x0, algorithm):
        self.iter += 1
        self._set_point(x0)
        if algorithm == 'phenomena':
            x1 = self._run_phenomena()
        elif algorithm == 'phenomena modular':
            x1 = self._run_phenomena_modular()
        elif algorithm == 'sequential modular':
            x1 = self._run_sequential()
        elif algorithm == 'inside out':
            x1 = self._run_inside_out()
        else:
            raise RuntimeError(f'invalid algorithm {algorithm!r}')
        x1 = self._new_point(x1)
        if self._convergence_analysis_mode: self._tracked_points[self.iter] = x1
        return x1
    
    # %% Inside-out simulation
    
    def _run_inside_out(self):
        T = np.zeros(self.N_stages)
        x = np.zeros((self.N_stages, self._N_chemicals))
        y = x.copy()
        gamma = x.copy()
        K = x.copy()
        dlogK_dTinv = x.copy()
        hV = T.copy()
        hL = T.copy()
        CV = T.copy()
        CL = T.copy()
        P = self.P
        f_gamma = self._gamma
        mixture = self._eq_thermo.mixture
        H = mixture.H
        C = mixture.Cn
        for i, stage in enumerate(self.stages):
            Pi = P[i]
            f = lambda Tinv: np.log(stage.K_model(stage.y, stage.x, 1 / Tinv[0], Pi))
            stage.partition._run_decoupled_KTvle(P=Pi)
            dlogK_dTinv[i] = approx_derivative(f, 1 / stage.T)[:, 0]
            T[i] = Ti = stage.T
            x[i] = xi = stage.x
            y[i] = yi = stage.y
            K[i] = stage.K
            gamma[i] = f_gamma(xi, Ti)
            hV[i] = H('g', yi, Ti, Pi)
            hL[i] = H('l', xi, Ti, Pi)
            CV[i] = C('g', yi, Ti, Pi)
            CL[i] = C('l', xi, Ti, Pi)
        V, L = MESH.bulk_vapor_and_liquid_flow_rates(
            hL, hV,
            self._neg_asplit, self._neg_bsplit, 
            self._top_split, self._bottom_split, 
            self.N_stages, self._feed_and_invariable_enthalpies, 
            self._total_feed_flows,
            self._specified_variables,
            self._specified_values,
            self._bulk_feed,
        )
        SC = SurrogateColumn(
            self.N_stages, 
            self._N_chemicals,
            T, x, y, 
            gamma, K, dlogK_dTinv,
            hV, hL,
            CV, CL,
            self._specified_variables,
            self._specified_values,
            self._neg_asplit,
            self._neg_bsplit,
            self._top_split,
            self._bottom_split,
            self._feed_and_invariable_enthalpies,
            self.feed_flows,
            self._total_feed_flows,
            self._bulk_feed,
            self._point_shape,
        )
        Kb = SC.Kb
        algorithm = self.inside_loop_algorithm
        Sb = Kb * V / L
        if algorithm == 'phenomena + simultaneous correction':
            f = SC.phenomena_iter
            Sb = flx.fixed_point(
                f, Sb, 
                xtol=self.tolerance, 
                rtol=self.relative_tolerance, 
                maxiter=self.maxiter,
                checkiter=False,
                checkconvergence=False,
            )
            f = SC.residuals
            jac = lambda logSb1: approx_derivative(f, logSb1)
            logSb1, *self._inside_info = fsolve(
                f, np.log(Sb + 1), fprime=jac, full_output=True, 
                maxfev=self.maxiter, xtol=self.relative_tolerance / 1000
            )
            Sb = np.exp(logSb1) - 1
        elif algorithm == 'phenomena':
            f = SC.phenomena_iter
            Sb = flx.fixed_point(
                f, Sb, 
                xtol=self.tolerance / 1000, 
                rtol=self.relative_tolerance / 1000, 
                maxiter=self.maxiter,
                checkiter=False,
                checkconvergence=False,
            )
        elif algorithm == 'simultaneous correction': 
            f = SC.residuals
            jac = lambda logSb1: approx_derivative(f, logSb1)
            logSb1, *self._inside_info = fsolve(
                f, np.log(Sb + 1), fprime=jac, full_output=True, 
                maxfev=self.maxiter, xtol=self.relative_tolerance / 1000,
            )
            Sb = np.exp(logSb1) - 1
        else:
            raise RuntimeError('invalid inside loop algorithm')
        return SC.Sb_to_point(Sb)
    
    # %% Normal simulation
    
    def _run_phenomena(self):
        if self._has_vle:
            decomp = self.vle_decomposition
            if decomp == 'bubble point':
                self.update_bubble_point()
                self.update_energy_balance_phase_ratios()
                for i in self.stages: i._update_separation_factors()
                separation_factors = np.array([i.S for i in self._S_stages])
                self.update_flow_rates(separation_factors, update_B=True)
            elif decomp == 'sum rates':
                self.update_pseudo_vle()
                separation_factors = np.array([i.S for i in self._S_stages])
                self.update_flow_rates(separation_factors, update_B=True)
                self.update_energy_balance_temperatures()
            elif decomp == 'bubble point Wang-Henke':
                mol_liq = self.get_liquid_flow_rates()
                xs = mol_liq / mol_liq.sum(axis=1, keepdims=True)
                Ks, ys, Ts = self.estimate_bubble_point(xs)
                Vs, Ls = self.estimate_bulk_vapor_and_liquid_flow_rates(xs, ys, Ts)
                Vs = Vs[:, None]
                Ls = Ls[:, None]
                xs = self.estimate_liquid_composition(Ks, Vs, Ls)
                self.update_WangHenke(ys * Vs, Ts, xs * Ls)
            else:
                raise NotImplementedError(f'{decomp!r} decomposition not implemented')
        elif self._has_lle:
            self.update_pseudo_lle()
            self.update_energy_balance_temperatures()
        else:
            raise RuntimeError('unknown equilibrium phenomena')
    
    def _run_phenomena_modular(self):
        self._run_sequential()
        stages = self.stages
        xs = np.array([i.x for i in stages])
        ys = np.array([i.y for i in stages])
        Ts = np.array([i.T for i in stages])
        
        # Energy balance
        Vs, Ls = self.estimate_bulk_vapor_and_liquid_flow_rates(xs, ys, Ts)
        Vs = Vs[:, None]
        Ls = Ls[:, None]
        Bs = Vs / Ls
        for i, stage in enumerate(self.stages):
            if stage.specified_variable != 'B': stage.B = Bs[i]
            vap, liq = stage.partition.outs
            if stage.specified_variable != 'T': stage.T = vap.T = liq.T = Ts[i]
        if getattr(self, 'tracking', False):
            self._collect_variables('energy')
        
        for i in self.stages: i._update_separation_factors()
        separation_factors = np.array([i.S for i in self._S_stages])
        self.update_flow_rates(separation_factors, update_B=True)
    
    def _run_sequential(self):
        *stages, last = self.stages
        for i in stages: i._run()
        last._run()
        for i in reversed(stages[1:]): i._run()
    
    def _simultaneous_correction(self, x, method):
        shape = x.shape
        if self._convergence_analysis_mode:
            self._tracked_algorithms.append(
                (self.iter + 1, 'simultaneous correction')
            )
            def f(x):
                x = x.reshape(shape)
                self.iter += 1
                try: self._tracked_points[self.iter] = x
                except: pass
                return self._residuals(x).flatten()
        else:
            f = lambda x: self._residuals(x.reshape(shape)).flatten()
        
        jac = lambda x: MESH.create_block_tridiagonal_matrix(*self._jacobian(x.reshape(shape)))
        try: 
            if method == 'trf':
                self._result = result = least_squares(f,
                    x0=x.flatten(),
                    jac=lambda x: MESH.create_block_tridiagonal_matrix(*self._jacobian(x.reshape(shape))),
                    bounds=(0, np.inf),
                    method='trf',
                    max_nfev=self.maxiter,
                    xtol=self.tolerance,
                    loss='cauchy',
                    x_scale='jac',
                    tr_solver='lsmr',
                    jac_sparsity=MESH.create_block_tridiagonal_matrix(*self._jacobian(x)),
                )
                x = result.x
            elif method == 'hybr':
                x, *self._simultaneous_correction_info = fsolve(
                    f, x.flatten(), fprime=jac, full_output=True, 
                    maxfev=self.maxiter, xtol=self.tolerance
                )
                
            else:
                raise ValueError(f'invalid simultaneous correction method {method!r}')
        except: pass
        else:
            x[x < 0] = 0
            x = x.reshape(shape)
            try: result = self._best_result
            except: self._best_result = result
            else:
                r = self._objective(x)
                if result.r < r: x = result.x
                else: self._best_result = IterationResult(None, x, r)
        return x
    
    
    # %% Initial guess
    
    def hot_start_collapsed_stages(self,
            all_stages, feed_stages, stage_specifications,
            top_side_draws, bottom_side_draws,
        ):
        raise NotImplementedError('not yet in BioSTEAM yet')
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
        self._eq_index = index = ms.chemicals.get_index(IDs)
        for i in stages: 
            i._eq_index = index
            i._N_chemicals = N_chemicals
        self.feed_flows = feed_flows = np.zeros([N_stages, N_chemicals])
        self.feed_enthalpies = feed_enthalpies = np.zeros(N_stages)
        for feed, stage in zip(feeds, feed_stages):
            feed_flows[stage, :] += feed.mol[index]
            feed_enthalpies[stage] += feed.H
        self._total_feed_flows = feed_flows.sum(axis=1)
        self._iter_args = (feed_flows, self._neg_asplit, self._neg_bsplit, self.N_stages)
        # feed_stages = [(i if i >= 0 else N_stages + i) for i in self.feed_stages]
        # stage_specifications = {(i if i >= 0 else N_stages + i): j for i, j in self.stage_specifications.items()}
        # top_side_draws = {(i if i >= 0 else N_stages + i): j for i, j in self.top_side_draws.items()}
        # bottom_side_draws = {(i if i >= 0 else N_stages + i): j for i, j in self.bottom_side_draws.items()}
        # self.key_stages = set([*feed_stages, *stage_specifications, *top_side_draws, *bottom_side_draws])
        self._bulk_feed = feed_flows.sum()
        # Reformulate flow rate specifications
        if (eq == 'vle'
            and stages[-1].specified_variable == 'F'
            and stages[0].specified_variable == 'B'):
            self._RF_spec = True
            last = stages[-1]
            last._bulk_feed = last.partition._bulk_feed = self._bulk_feed
            self._top_bulk_feed = feed_flows[0].sum()
        else:
            self._RF_spec = False
        N_chemicals = len(index)
        if top_chemicals:
            top_side_draws = self.top_side_draws
            n = len(top_chemicals)
            b = np.ones([N_stages, n])
            c = self._neg_asplit[1:]
            d = np.zeros([N_stages, n])
            for feed, stage in zip(feeds, feed_stages):
                d[stage] += feed.imol[top_chemicals]
            top_flow_rates = MESH.solve_right_bidiagonal_matrix(b, c, d)
            for partition, flows in zip(partitions, top_flow_rates):
                partition.outs[0].imol[top_chemicals] = flows
        if bottom_chemicals:
            bottom_side_draws = self.bottom_side_draws
            a = self._neg_bsplit[:-1]
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
        other_chemicals = top_chemicals + bottom_chemicals
        self._noneq_index = self.thermo.chemicals.indices(other_chemicals)
        self._noneq_thermo = self.thermo.subset(other_chemicals)
        self._eq_thermo = self.thermo.subset(IDs)
        Hother = self._noneq_thermo.mixture.H
        P = self.P
        N_stages = self.N_stages
        variables = ''
        self._specified_values = values = np.zeros(N_stages)
        invariable_enthalpies = np.zeros(N_stages)
        other_index = self._noneq_index
        for n, stage in enumerate(self.stages):
            Pi = P[n]
            variable = stage.specified_variable
            if variable == 'T':
                variables += 'B'
                values[n] = stage.B
            else:
                variables += variable
                values[n] = getattr(stage, variable)
            partition = stage.partition
            vap, liq = partition.outs
            H_in = sum([
                Hother(i.phase, i.mol[other_index], i.T, Pi)
                + Hother(i.phase, i.mol[other_index], i.T, Pi)
                for i in stage.ins
            ])
            H_out = sum([
                Hother(i.phase, i.mol[other_index], i.T, Pi)
                + Hother(i.phase, i.mol[other_index], i.T, Pi)
                for i in partition.outs
            ])
            invariable_enthalpies[n] += H_out - H_in
        self._feed_and_invariable_enthalpies = invariable_enthalpies + feed_enthalpies
        self._specified_variables = variables
        if (self.use_cache 
            and all([i.IDs == IDs for i in partitions])): # Use last set of data
            pass
        else:
            for i in partitions: i.IDs = IDs
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
                for i in self.stages: i._update_separation_factors()
                self.update_flow_rates(np.array([i.S for i in self._S_stages]), update_B=False)
            elif eq == 'lle':
                lle = ms.lle
                T = ms.T
                lle(T, top_chemical=top_chemical)
                K = lle._K
                phi = lle._phi
                B = inf if phi == 1 else phi / (1 - phi)
                y_mol = ms.imol['L', IDs]
                x_mol = ms.imol['l', IDs]
                y = y_mol / y_mol.sum()
                f_gamma = self.thermo.Gamma([chemicals[i] for i in IDs])
                gamma_y = f_gamma(y, T)
                for i in partitions: 
                    i.B = B
                    i.T = T
                    i.K = K
                    i.gamma_y = gamma_y
                    for j in i.outs: j.T = T
                top_flows = np.ones((N_stages, N_chemicals)) * y_mol
                bottom_flows = np.ones((N_stages, N_chemicals)) * x_mol
                self.set_all_flow_rates(top_flows, bottom_flows)
            elif self.stage_specifications:
                P = self.P[N_stages // 2]
                dp = ms.dew_point_at_P(P=P, IDs=IDs)
                T_bot = dp.T
                bp = ms.bubble_point_at_P(P=P, IDs=IDs)
                T_top = bp.T
                dT_stage = (T_bot - T_top) / N_stages
                Ts = np.array([T_top + i * dT_stage for i in range(N_stages)])
                z = bp.z.copy()
                z[z == 0] = 1e-32
                x = dp.x.copy()
                x[x == 0] = 1e-32
                K_dew = dp.z / dp.x
                K_bubble = bp.y / bp.z
                K_dew[K_dew > 1e9] = 1e9
                K_bubble[K_bubble > 1e9] = 1e9
                K_dew[K_dew < 1e-9] = 1e-9
                K_bubble[K_bubble < 1e-9] = 1e-9
                dK_stage = (K_bubble - K_dew) / N_stages
                top_flows = np.ones((N_stages, N_chemicals))
                bottom_flows = np.ones((N_stages, N_chemicals))
                for i, partition in enumerate(partitions):
                    partition.T = T = T_top + i * dT_stage
                    partition.K = K_dew + i * dK_stage
                    a = i / N_stages
                    b = 1 - a
                    x = a * bp.x + b * dp.x
                    x /= x.sum()
                    y = a * bp.y + b * dp.y
                    y /= y.sum()
                    bottom_flows[i] = partition.x = x
                    top_flows[i] = partition.y = y
                    for s in partition.outs: s.T = T
                xs = np.array([i.x for i in partitions])
                ys = np.array([i.y for i in partitions])
                Vs, Ls = self.estimate_bulk_vapor_and_liquid_flow_rates(xs, ys, Ts)
                phase_ratios = Vs / Ls
                for partition, B in zip(partitions, phase_ratios):
                    if partition.specified_variable != 'B': partition.B = B
                top_flows *= Vs[:, None]
                bottom_flows *= Ls[:, None]
                self.set_all_flow_rates(top_flows, bottom_flows)
            else:
                vle = ms.vle
                P = self.P[N_stages // 2]
                vle(H=ms.H, P=P)
                L_mol = ms.imol['l', IDs]
                x = L_mol / L_mol.sum()
                V_mol = ms.imol['g', IDs]
                y = V_mol / V_mol.sum()
                K = y / x
                phi = ms.V
                B = phi / (1 - phi)
                T = ms.T
                for P, partition in zip(self.P, partitions):
                    partition.T = T
                    partition.B = B
                    for i in partition.outs: i.T = T
                    partition.K = K
                    partition.fgas = P * y
                    partition.y = y
                    partition.x = x
                    for s in partition.outs: s.empty()
                top_flows = np.ones((N_stages, N_chemicals)) * V_mol
                bottom_flows = np.ones((N_stages, N_chemicals)) * L_mol
                self.set_all_flow_rates(top_flows, bottom_flows)
        if eq == 'vle' and self.vle_decomposition is None: self.default_vle_decomposition()
        self._gamma = self._eq_thermo.Gamma(self._eq_thermo.chemicals)
        self._H_magnitude = 10 * sum([i.mixture.Cn('l', i.mol, i.T, i.P) for i in self.ins])
        self.attempt = 0
        self._mean_residual = np.inf
        self._best_result = empty = IterationResult(None, None, np.inf)
        self._point_shape = (N_stages, 2 * N_chemicals + 1)
        record = self.iteration_memory * [empty]
        x = self._get_point()
        record[0] = IterationResult(1, x, self._objective(x))
        self._iteration_record = record = deque(record)
        return x
    
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
                        self._asplit[index], 
                        self._bsplit[index],
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
                Cv, Cl, Hv, Hl, self._asplit, self._bsplit,
                N_stages, self.feed_enthalpies
            )
            if not np.isfinite(dTs).all():
                breakpoint()
                dTs = MESH.temperature_departures(
                    Cv, Cl, Hv, Hl, self._asplit, self._bsplit,
                    N_stages, self.feed_enthalpies
                )
        return dTs
    
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
    
    def update_energy_balance_phase_ratios(self):
        stages = self.stages
        xs = np.array([i.x for i in stages])
        ys = np.array([i.y for i in stages])
        Ts = np.array([i.T for i in stages])
        
        # Energy balance
        Vs, Ls = self.estimate_bulk_vapor_and_liquid_flow_rates(xs, ys, Ts)
        Vs = Vs[:, None]
        Ls = Ls[:, None]
        Bs = Vs / Ls
        for stage, B in zip(stages, Bs):
            if stage.specified_variable != 'B': stage.B = B
        if getattr(self, 'tracking', False):
            self._collect_variables('energy')
    
    def estimate_bulk_vapor_and_liquid_flow_rates(self, xs, ys, Ts):
        Hvle = self._eq_thermo.mixture.H
        N_stages = self.N_stages
        return MESH.bulk_vapor_and_liquid_flow_rates(
                np.array([Hvle('l', i.x, i.T, i.P) for i in self.stages]), 
                np.array([Hvle('g', i.y, i.T, i.P) for i in self.stages]), 
                self._neg_asplit, self._neg_bsplit, 
                self._top_split, self._bottom_split, 
                N_stages, self._feed_and_invariable_enthalpies, 
                self._total_feed_flows,
                self._specified_variables,
                self._specified_values,
                self._bulk_feed,
            )
    
    # %% Material balance convergence
    
    def _feed_flows_and_conversion(self):
        feed_flows = self.feed_flows.copy()
        partition = self.partitions
        index = self._eq_index
        for i in self.stage_reactions: 
            p = partition[i]
            dmol = p.dmol
            for n, j in enumerate(index): feed_flows[i, n] += dmol[j]
        return feed_flows
    
    def set_all_flow_rates(self, top_flows, bottom_flows):
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        index = self._eq_index
        RF_spec = self._RF_spec
        for i in range_stages:
            stage = stages[i]
            partition = stage.partition
            s_top, s_bot = partition.outs
            t = top_flows[i]
            mask = t < 0
            bulk_t = t.sum()
            if mask.any():
                t[mask] = 0
                tsum = t.sum()
                if tsum: t *= bulk_t / tsum
            b = bottom_flows[i]
            mask = b < 0
            bulk_b = b.sum()
            if bulk_b < 0:
                bottom_flows[i] = bottom_flows[i-1]
                bulk_b = b.sum()
            elif mask.any():
                b[mask] = 0
                bsum = b.sum()
                if bsum != 0: b *= bulk_b / bsum
            for i in stage.splitters: i._run()
            if stage.specified_variable == 'B':
                if stage.B == 0:
                    t[:] = 0
                else:
                    t *= stage.B * bulk_b / bulk_t
            elif RF_spec and i == 1:
                first = stages[0]
                second = stage
                last = stages[-1]
                F_feed = self._bulk_feed
                F_bot = last.F * F_feed
                F_dist = F_feed - F_bot
                F_reflux = F_dist / first.B
                F_second_vap = (F_reflux + F_dist - self._top_bulk_feed) / (1 - second.top_split)
                F_second_liq = bulk_b
                second.B = F_second_vap / F_second_liq
                bottom_flows[-1] *= F_bot / bottom_flows[-1].sum()
                t *= F_second_vap / bulk_t
            s_top.mol[index] = t
            s_bot.mol[index] = b
    
    def set_flow_rates(self, bottom_flows, update_B=True):
        stages = self.stages
        N_stages = self.N_stages
        range_stages = range(N_stages)
        feed_flows = self.feed_flows
        index = self._eq_index
        if self.stage_reactions:
            feed_flows = self._feed_flows_and_conversion()
        top_flows = MESH.top_flows_mass_balance(
            bottom_flows, feed_flows, self._asplit, self._bsplit, 
            self.N_stages
        )
        f = PhasePartition.F_relaxation_factor
        if f != 0: raise NotImplementedError('F relaxation factor')
        RF_spec = self._RF_spec
        for i in range_stages:
            stage = stages[i]
            partition = stage.partition
            s_top, s_bot = partition.outs
            t = top_flows[i]
            mask = t < 0
            bulk_t = t.sum()
            if mask.any():
                t[mask] = 0
                tsum = t.sum()
                if tsum: t *= bulk_t / tsum
            b = bottom_flows[i]
            mask = b < 0
            bulk_b = b.sum()
            if bulk_b < 0:
                bottom_flows[i] = bottom_flows[i-1]
                bulk_b = b.sum()
            elif mask.any():
                b[mask] = 0
                bsum = b.sum()
                if bsum != 0: b *= bulk_b / bsum
            for i in stage.splitters: i._run()
            if stage.specified_variable == 'B':
                if stage.B == 0:
                    t[:] = 0
                else:
                    t *= stage.B * bulk_b / bulk_t
            elif RF_spec and i == 1:
                first = stages[0]
                second = stage
                last = stages[-1]
                F_feed = self._bulk_feed
                F_bot = last.F * F_feed
                F_dist = F_feed - F_bot
                F_reflux = F_dist / first.B
                F_second_vap = (F_reflux + F_dist - self._top_bulk_feed) / (1 - second.top_split)
                F_second_liq = bulk_b
                second.B = F_second_vap / F_second_liq
                bottom_flows[-1] *= F_bot / bottom_flows[-1].sum()
                t *= F_second_vap / bulk_t
            elif update_B:
                stage.B = bulk_t / bulk_b
            s_top.mol[index] = t
            s_bot.mol[index] = b
    
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
            denominator = S_old - S
            denominator[denominator == 0] = 1
            x = S_old / denominator
            x = 0.1 * (x[(x < 1) & (x > 0)]).min()
            separation_factors = S * x + (1 - x) * S_old
        for stage, S in zip(self._S_stages, separation_factors): stage.S = S
        flows = self.run_mass_balance()
        if not np.isfinite(flows).all():
            raise RuntimeError('infeasible values in flow rates')
        self.set_flow_rates(flows, update_B)
        if getattr(self, 'tracking', False):
            self._collect_variables('material')
    
    # %% Phenomena convergence
    
    def update_bubble_point(self):
        stages = self.stages
        P = self.P
        if self.stage_reactions:
            self.update_liquid_holdup() # Finds liquid volume at each stage
            for n, stage in enumerate(stages):
                partition = stage.partition
                partition._run_decoupled_KTvle(P=P[n])
                T = partition.T
                for i in (partition.outs + stage.outs): i.T = T
                if partition.reaction: 
                    partition._run_decoupled_reaction(P=P)
        else:
            for n, stage in enumerate(stages):
                partition = stage.partition
                partition._run_decoupled_KTvle(P=P[n])
                T = partition.T
                for i in (partition.outs + stage.outs): i.T = T
        if getattr(self, 'tracking', False):
            self._collect_variables('vle_phenomena')
    
    def update_pseudo_lle(self):
        stages = self.stages
        if 'K' in self.partition_data:
            for i in self.stages:
                i.K = i.partition.partition_data['K']
            self.update_mass_balance()
        else:
            stages = self._S_stages
            def psuedo_equilibrium(flow_rates):
                self.set_flow_rates(flow_rates, update_B=False)
                for n, i in enumerate(stages): 
                    i.partition._run_decoupled_Kgamma()
                    i._update_separation_factors()
                return self.run_mass_balance()
            
            separation_factors = flx.fixed_point(
                psuedo_equilibrium, self.run_mass_balance(), 
                xtol=self.tolerance,
                rtol=self.relative_tolerance,
                checkiter=False,
                checkconvergence=False,
            )
        for i in stages: 
            mixer = i.mixer
            partition = i.partition
            mixer.outs[0].mix_from(mixer.ins, energy_balance=False)
            partition._run_decoupled_B()
            i._update_separation_factors()
        separation_factors = np.array([i.S for i in stages])
        self.update_flow_rates(separation_factors, update_B=False)
        if getattr(self, 'tracking', False):
            self._collect_variables('lle_phenomena')
    
    def update_pseudo_vle(self, separation_factors):
        P = self.P
        if 'K' not in self.partition_data:
            for n, i in enumerate(self.stages): 
                i.partition._run_decoupled_Kfgas(P=P[n])
                i._update_separation_factors()
        if getattr(self, 'tracking', False):
            raise NotImplementedError('tracking sum rates decomposition not implemented')
        
# %% Wang-Henke

    def get_liquid_flow_rates(self):
        index = self._eq_index
        return np.array([i.outs[1].mol[index] for i in self.partitions])
    
    def get_bulk_vapor_flow_rates(self):
        index = self._eq_index
        return np.array([i.outs[0].mol[index].sum() for i in self.partitions])
    
    def estimate_bubble_point(self, xs):
        bubble_point = tmo.equilibrium.BubblePoint(thermo=self._eq_thermo)
        P = self.P
        Tys = [bubble_point.solve_Ty(x, P) for x in xs]
        Ts = np.zeros(self.N_stages)
        Ks = xs.copy()
        ys = xs.copy()
        for i, (T, y) in enumerate(Tys):
            Ks[i] = y / xs[i]
            ys[i] = y
            Ts[i] = T
        return Ks, ys, Ts
    
    def estimate_liquid_composition(self, Ks, Vs, Ls):
        xs = MESH.liquid_compositions(Vs, Ls, Ks, self.feed_flows, self._neg_asplit, self._neg_bsplit, self.N_stages)
        xs[xs < 0] = 1e-16
        return xs / xs.sum(axis=1, keepdims=True)
    
    def update_WangHenke(self, mol_vap, Ts, mol_liq):
        x = np.zeros(self._point_shape)
        N_chemicals = self._N_chemicals
        x[:, :N_chemicals] = mol_vap
        x[:, N_chemicals] = Ts
        x[:, -N_chemicals:] = mol_liq
        self._set_point(x)
        
# %% Convergence analysis

    def convergence_analysis(
            self, 
            iterations=None, 
            algorithm=None, 
            x0=None, xf=None,
            fillsteps=1,
            yticks=None,
            xticks=None,
            method=None,
            **kwargs,
        ):
        x0 = self.hot_start() if x0 is None else x0
        if xf is None:
            self._run()
            xf = self._get_point()
            self._set_point(x0)
        diff = (xf - x0)
        maxdistance = np.sqrt((diff * diff).sum())
        if self.vle_decomposition is None: self.default_vle_decomposition()
        if iterations is None: iterations = self.maxiter
        self._convergence_analysis_mode = True
        try:
            self._tracked_points = points = np.zeros([iterations + 1, self.N_stages, self._N_chemicals * 2 + 1])
            self._tracked_algorithms = algorithms = []
            self._get_point(points[0])
            if algorithm is None:
                self._run()
                iterations = min(self.iter, iterations)
            elif algorithm == 'simultaneous correction':
                shape = x0.shape
                self.iter = 0
                def f(x):
                    x = x.reshape(shape)
                    self.iter += 1
                    points[self.iter] = x
                    return self._residuals(x).flatten()
                
                if method == 'trf':
                    least_squares(f,
                        x0=x0.flatten(),
                        jac=lambda x: MESH.create_block_tridiagonal_matrix(*self._jacobian(x.reshape(shape))),
                        bounds=(0, np.inf),
                        method='trf',
                        max_nfev=iterations,
                        x_scale='jac',
                        tr_solver='lsmr',
                        jac_sparsity=MESH.create_block_tridiagonal_matrix(*self._jacobian(x0)),
                    )
                elif method == 'hybr':
                    jac = lambda x: MESH.create_block_tridiagonal_matrix(*self._jacobian(x.reshape(shape)))
                    try: fsolve(f, x0.flatten(), fprime=jac, maxfev=iterations)
                    except: pass
                else:
                    raise ValueError('invalid method')
                iterations = min(self.iter, iterations)
            else:
                self.iter = 0
                def f(x):
                    self.iter += 1
                    self._set_point(x)
                    if algorithm == 'phenomena':
                        self._run_phenomena()
                    elif algorithm == 'sequential modular':
                        self._run_sequential()
                    elif algorithm == 'phenomena modular':
                        self._run_phenomena_modular()
                    elif algorithm == 'inside out':
                        try:
                            points[self.iter] = self._run_inside_out()
                        except:
                            pass
                        return points[self.iter]
                    else:
                        raise ValueError('invalid algorithm')
                    return self._get_point(points[self.iter])
                
                if method == 'wegstein':
                    solver = flx.wegstein
                elif method is None or method == 'fixed-point':
                    solver = flx.fixed_point
                else:
                    raise ValueError('unknown method')
                try:
                    solver(f, x0, xtol=0, 
                           maxiter=iterations-1, 
                           checkconvergence=False, 
                           checkiter=False,
                           **kwargs)
                except:
                    pass
        finally:
            self._convergence_analysis_mode = False
        iterations -= 1
        corrections = np.diff(points, axis=0)
        stepsize = 1 / fillsteps
        shape = (iterations, fillsteps)
        distances = np.zeros(shape)
        residuals = np.zeros(shape)
        iteration = np.zeros(shape)
        f = self._objective
        for i in range(iterations):
            x0 = points[i]
            dx = corrections[i]
            for j in range(fillsteps):
                t = j * stepsize
                x = x0 + t * dx
                diff = xf - x
                iteration[i, j] = i + t
                distances[i, j] = np.sqrt((diff * diff).sum())
                residuals[i, j] = f(x)
        colors = bst.utils.GG_colors
        blue = colors.blue.RGBn
        red = colors.red.RGBn
        cmap = clr.LinearSegmentedColormap.from_list(
            'blue2red',
            [blue, red],
            N=256
        )
        mindistance = 1e-1
        distances[distances > maxdistance] = maxdistance
        distances[distances < mindistance] = mindistance
        distances = np.log(distances) 
        distances -= np.log(mindistance)
        distances /= np.log(maxdistance) - np.log(mindistance)
        fig = plt.figure()
        ax = plt.gca()
        plt.scatter(
            iteration.flatten(),
            np.log(residuals.flatten()),
            c=cmap(distances.flatten()),
        )
        if yticks is not None: 
            plt.yticks(yticks)
            lb, ub = yticks[0], yticks[-1]
            plt.ylim([lb, ub])
        else:
            lb, ub = plt.ylim()
        if xticks is not None:
            plt.xticks(xticks)
            plt.xlim([xticks[0], xticks[-1]])
        sm = cm.ScalarMappable(cmap=cmap)    
        cbar = fig.colorbar(sm, ax=ax, ticks=[0, 1])
        
        yloc = lb + (ub - lb) * 1.05
        if algorithms:
            shorthand = {
                'phenomena': 'P',
                'phenomena modular': 'PM',
                'sequential modular': 'SM',
                'simultaneous correction': 'SC',
                'inside out': 'IO',
            }
            for iteration, algorithm in algorithms:
                plt.axvline(iteration, color='silver', zorder=-1)
                ax.text(
                    iteration, yloc, shorthand[algorithm], 
                    ha='center',
                    fontdict=dict(
                        fontname='arial', 
                        size=10
                    )
                )
        else:
            ax.text(
                np.mean(plt.xlim()), yloc, algorithm, 
                ha='center',
                fontdict=dict(
                    fontname='arial', 
                    size=10
                )
            )
        plt.xlabel('Iteration')
        plt.ylabel('log residual')
        cbar.ax.set_ylabel('distance from steady state')
        
        # cbar_ax.set_title()
            
        
        
        
    