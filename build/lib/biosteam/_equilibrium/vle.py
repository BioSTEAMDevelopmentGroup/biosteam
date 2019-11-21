# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:40:05 2019

@author: yoelr
"""
from flexsolve import bounded_wegstein, wegstein, aitken, \
                      bounded_aitken, IQ_interpolation
from numba import njit
import numpy as np

__all__ = ('VLE', 'V_2N', 'V_3N', 'V_error')

@njit
def V_error(V, zs, Ks):
    """Vapor fraction error."""
    return (zs*(Ks-1.)/(1.+V*(Ks-1.))).sum()

@njit
def V_2N(zs, Ks):
    """Solution for 2 component flash vessel."""
    z1, z2 = zs
    K1, K2 = Ks
    return (-K1*z1 - K2*z2 + z1 + z2)/(K1*K2*z1 + K1*K2 *
                                       z2 - K1*z1 - K1*z2
                                       - K2*z1 - K2*z2 + z1 + z2)
@njit    
def V_3N(zs, Ks):
    """Solution for 3 component flash vessel."""
    z1, z2, z3 = zs
    K1, K2, K3 = Ks
    return (-K1*K2*z1/2 - K1*K2*z2/2 - K1*K3*z1/2 - K1*K3*z3/2 + K1*z1 + K1*z2/2 + K1*z3/2 - K2*K3*z2/2 - K2*K3*z3/2 + K2*z1/2 + K2*z2 + K2*z3/2 + K3*z1/2 + K3*z2/2 + K3*z3 - z1 - z2 - z3 - (K1**2*K2**2*z1**2 + 2*K1**2*K2**2*z1*z2 + K1**2*K2**2*z2**2 - 2*K1**2*K2*K3*z1**2 - 2*K1**2*K2*K3*z1*z2 - 2*K1**2*K2*K3*z1*z3 + 2*K1**2*K2*K3*z2*z3 - 2*K1**2*K2*z1*z2 + 2*K1**2*K2*z1*z3 - 2*K1**2*K2*z2**2 - 2*K1**2*K2*z2*z3 + K1**2*K3**2*z1**2 + 2*K1**2*K3**2*z1*z3 + K1**2*K3**2*z3**2 + 2*K1**2*K3*z1*z2 - 2*K1**2*K3*z1*z3 - 2*K1**2*K3*z2*z3 - 2*K1**2*K3*z3**2 + K1**2*z2**2 + 2*K1**2*z2*z3 + K1**2*z3**2 - 2*K1*K2**2*K3*z1*z2 + 2*K1*K2**2*K3*z1*z3 - 2*K1*K2**2*K3*z2**2 - 2*K1*K2**2*K3*z2*z3 - 2*K1*K2**2*z1**2 - 2*K1*K2**2*z1*z2 - 2*K1*K2**2*z1*z3 + 2*K1*K2**2*z2*z3 + 2*K1*K2*K3**2*z1*z2 - 2*K1*K2*K3**2*z1*z3 - 2*K1*K2*K3**2*z2*z3 - 2*K1*K2*K3**2*z3**2 + 4*K1*K2*K3*z1**2 + 4*K1*K2*K3*z1*z2 + 4*K1*K2*K3*z1*z3 + 4*K1*K2*K3*z2**2 + 4*K1*K2*K3*z2*z3 + 4*K1*K2*K3*z3**2 + 2*K1*K2*z1*z2 - 2*K1*K2*z1*z3 - 2*K1*K2*z2*z3 - 2*K1*K2*z3**2 - 2*K1*K3**2*z1**2 - 2*K1*K3**2*z1*z2 - 2*K1*K3**2*z1*z3 + 2*K1*K3**2*z2*z3 - 2*K1*K3*z1*z2 + 2*K1*K3*z1*z3 - 2*K1*K3*z2**2 - 2*K1*K3*z2*z3 + K2**2*K3**2*z2**2 + 2*K2**2*K3**2*z2*z3 + K2**2*K3**2*z3**2 + 2*K2**2*K3*z1*z2 - 2*K2**2*K3*z1*z3 - 2*K2**2*K3*z2*z3 - 2*K2**2*K3*z3**2 + K2**2*z1**2 + 2*K2**2*z1*z3 + K2**2*z3**2 - 2*K2*K3**2*z1*z2 + 2*K2*K3**2*z1*z3 - 2*K2*K3**2*z2**2 - 2*K2*K3**2*z2*z3 - 2*K2*K3*z1**2 - 2*K2*K3*z1*z2 - 2*K2*K3*z1*z3 + 2*K2*K3*z2*z3 + K3**2*z1**2 + 2*K3**2*z1*z2 + K3**2*z2**2)**0.5/2)/(K1*K2*K3*z1 + K1*K2*K3*z2 + K1*K2*K3*z3 - K1*K2*z1 - K1*K2*z2 - K1*K2*z3 - K1*K3*z1 - K1*K3*z2 - K1*K3*z3 + K1*z1 + K1*z2 + K1*z3 - K2*K3*z1 - K2*K3*z2 - K2*K3*z3 + K2*z1 + K2*z2 + K2*z3 + K3*z1 + K3*z2 + K3*z3 - z1 - z2 - z3)

class VLE:
    """Create a VLE object for solving VLE."""
    __slots__ = ('T', 'P', 'Q', 'V', '_stream', '_gamma',
                 '_dew_point', '_bubble_point', '_v', '_liquid_mol',
                 '_vapor_mol', '_index', '_massnet', '_compound',
                 '_update_V', '_mol', '_molnet', '_N', '_solve_V',
                 '_zs', '_Ks', '_Psat_gama', '_Psat_P')
    
    solver = staticmethod(IQ_interpolation)
    itersolver = staticmethod(aitken)
    T_tol = 0.00001
    P_tol = 0.1
    Q_tol = 0.1
    V_tol = 0.00001
    
    def __init__(self, stream):
        self.T = self.P = self.Q = self.V = 0
        self._stream = stream
        self._dew_point = stream._dew_point
        self._bubble_point = stream._bubble_point
        self._gamma = stream._gamma
        self._liquid_mol = stream.liquid_mol
        self._vapor_mol = stream.vapor_mol
    
    def __call__(self, species_IDs=None, LNK=None, HNK=None,
                 P=None, Q=None, T=None, V=None, x=None, y=None):
        """Perform vapor-liquid equilibrium.

        Parameters
        ----------
        Specify two:
            * **P:** Operating pressure (Pa)
            * **Q:** Energy input (kJ/hr)
            * **T:** Operating temperature (K)
            * **V:** Molar vapor fraction
            * **x:** Molar composition of liquid (for binary mixture)
            * **y:** Molar composition of vapor (for binary mixture)
        species_IDs = None : tuple, optional
                             IDs of species in equilibrium.
        LNK = None : tuple[str], optional
              Light non-keys that remain as a vapor (disregards equilibrium).
        LNK = None : tuple[str], optional
              Heavy non-keys that remain as a liquid (disregards equilibrium).

        """
        ### Decide what kind of equilibrium to run ###
        T_spec = T is not None
        P_spec = P is not None
        V_spec = V is not None
        Q_spec = Q is not None
        x_spec = x is not None
        y_spec = y is not None
        N_specs = (P_spec + Q_spec + T_spec + V_spec + x_spec + y_spec)
        assert N_specs == 2, ("must pass two and only two of the following specifications: "
                              "T, P, V, Q, x, y")
        self.setup(species_IDs, LNK, HNK)
        if T_spec:
            if P_spec:
                return self.TP(T, P)
            elif V_spec:
                return self.TV(T, V)
            elif Q_spec:
                return self.TQ(T, Q)
            elif x_spec:
                return self.Tx(T, np.asarray(x))
            else: # y_spec
                return self.Ty(T, np.asarray(y))
        elif P_spec:
            if V_spec:
                return self.PV(P, V)
            elif Q_spec:
                return self.PQ(P, Q)
            elif x_spec:
                return self.Px(P, np.asarray(x))
            else: # y_spec
                return self.Py(P, np.asarray(y))
        elif Q_spec:
            if y_spec:
                raise NotImplementedError('specification Q and y not implemented')
            elif x_spec:
                raise NotImplementedError('specification Q and x not implemented')
            else: # V_spec
                raise NotImplementedError('specification V and Q not implemented yet')
        else: # x_spec and y_spec
            raise ValueError("can only pass either 'x' or 'y' arguments, not both")
    
    def setup(self, species_IDs, LNK, HNK):
        liquid_mol = self._liquid_mol 
        vapor_mol = self._vapor_mol

        # Get flow rates
        mol = liquid_mol + vapor_mol
        notzero = mol > 0

        # Reused attributes
        species = self._stream._species
        speciesdct = species.__dict__
        indexdct = species._indexdct

        # Set up indices for both equilibrium and non-equilibrium species
        if species_IDs:
            index = [indexdct[specie] for specie in species_IDs]
            self._gamma.species = [speciesdct[ID] for ID in species_IDs]
        else:
            index = species._equilibrium_indices(notzero)
            cmps = species._compounds
            self._gamma.species = [cmps[i] for i in index]
        if LNK:
            LNK_index = [indexdct[i] for i in LNK]
        else:
            LNK_index = species._light_indices(notzero)
        if HNK:
            HNK_index = [indexdct[i] for i in HNK]
        else:
            HNK_index = species._heavy_indices(notzero)
        self._mol = mol[index]
        self._index = index

        # Set light and heavy keys
        vapor_mol[HNK_index] = 0
        vapor_mol[LNK_index] = mol[LNK_index]
        liquid_mol[LNK_index] = 0
        liquid_mol[HNK_index] = mol[HNK_index]
        
        self._N = N = len(index)
        if N == 1:
            self._compound = self._gamma._species[0]
            return 
        elif N == 2:
            self._solve_V = self._solve_V_2
        elif N == 3:
            self._solve_V = self._solve_V_3
        else:
            self._solve_V = self._solve_V_N
        self._massnet = self._stream.massnet
        
        # Get overall composition
        self._molnet = molnet = self._mol.sum()
        assert molnet != 0, 'empty stream cannot perform equilibrium'
        self._zs = self._mol/molnet

    ### Single component equilibrium case ###
        
    def _TP_compound(self, T, P):
        # Either liquid or gas
        if P < self._compound.VaporPressure(T):
            self._liquid_mol[self.index] = 0
            self._vapor_mol[self.index] = self._mol
        else:
            self._liquid_mol[self._index] = self._mol
            self._vapor_mol[self._index] = 0
    
    def _TV_compound(self, T, V):
        # Set vapor fraction
        self._stream.P = self._compound.VaporPressure(T)
        self._vapor_mol[self._index] = self._mol*V
        self._liquid_mol[self._index] = self._mol - self._vapor_mol[self._index]
        
    def _PV_compound(self, P, V):
        # Set vapor fraction
        self.T = self._compound.Tsat(P)
        self._vapor_mol[self._index] = self._mol*V
        self._liquid_mol[self._index] = self._mol - self._vapor_mol[self._index]
        
    def _PQ_compound(self, P, Q): 
        mol = self._mol
        index = self._index
        stream = self._stream
        vapor_mol = self._vapor_mol
        liquid_mol = self._liquid_mol
        Hin = Q + stream.H
        
        # Set temperature in equilibrium
        stream.T = self._compound.Tsat(P)
        
        # Check if super heated vapor
        vapor_mol[index] = mol
        liquid_mol[index] = 0
        H_dew = stream.H
        if Hin >= H_dew:
            stream.H = Hin
            return

        # Check if subcooled liquid
        vapor_mol[index] = 0
        liquid_mol[index] = mol
        H_bubble = stream.H
        if Hin <= H_bubble:
            stream.H = Hin
            return
        
        # Adjust vapor fraction accordingly
        V = (Hin - H_bubble)/(H_dew - H_bubble)
        vapor_mol[index] = mol*V
        liquid_mol[index] = mol - vapor_mol[index]
        
    def _TQ_compound(self, T, Q):
        index = self._index
        mol = self._mol
        stream = self._stream
        vapor_mol = self._vapor_mol
        liquid_mol = self._liquid_mol
        Hin = Q + stream.H
        
        # Set Pressure in equilibrium
        stream.P = self._compound.VaporPressure(T)
        
        # Check if super heated vapor
        vapor_mol[index] = mol
        liquid_mol[index] = 0
        H_dew = stream.H
        if Hin >= H_dew:
            stream.H = Hin
            return

        # Check if subcooled liquid
        vapor_mol[index] = 0
        liquid_mol[index] = mol
        H_bubble = stream.H
        if Hin <= H_bubble:
            stream.H = Hin
            return
        
        # Adjust vapor fraction accordingly
        V = (Hin - H_bubble)/(H_dew - H_bubble)
        vapor_mol[index] = mol*V
        liquid_mol[index] = mol - vapor_mol[index]
        
    def _lever_rule(self, x, y):
        split_frac = (self._zs[0]-x[0])/(y[0]-x[0])
        assert -0.00001 < split_frac < 1.00001, 'desired composition is infeasible'
        if split_frac > 1:
            split_frac = 1
        elif split_frac < 0:
            split_frac = 0
        self._vapor_mol[self._index] = self._molnet * split_frac * y
        self._liquid_mol[self._index] = self._mol - self._vapor_mol[self._index]
    
    def Tx(self, T, x):
        assert self._N == 2, 'number of species in equilibrium must be 2 to specify x'
        self._stream.P, y = self._bubble_point.solve_Px(x, T)
        self._lever_rule(x, y)
    
    def Px(self, P, x):
        assert self._N == 2, 'number of species in equilibrium must be 2 to specify x'
        self._stream.T, y = self._bubble_point.solve_Ty(x, P) 
        self._lever_rule(x, y)
        
    def Ty(self, T, y):
        assert self._N == 2, 'number of species in equilibrium must be 2 to specify y'
        self._stream.P, x = self._dew_point.solve_Px(y, T)
        self._lever_rule(x, y)
    
    def Py(self, P, y):
        assert self._N == 2, 'number of species in equilibrium must be 2 to specify y'
        self._stream.T, x = self._dew_point.solve_Ty(y, P) 
        self._lever_rule(x, y)
        
    def TP(self, T, P):
        stream = self._stream
        self.T = stream.T = T
        self.P = self._stream.P = P
        if self._N == 1: return self._TP_compound(T, P)
        # Setup bounderies
        P_dew, x_dew = self._dew_point.solve_Px(self._zs, T)
        P_bubble, y_bubble = self._bubble_point.solve_Py(self._zs, T)
        
        # Check if there is equilibrium
        if T >= P_dew:
            self._vapor_mol[self._index] = self._mol
            self._liquid_mol[self._index] = 0
        elif T <= P_bubble:
            self._vapor_mol[self._index] = 0
            self._liquid_mol[self._index] = self._mol
        else:
            # Guess composition in the vapor is a
            # weighted average of bubble/dew points
            self.V = V = self.V or (T - P_dew)/(P_bubble - P_dew)
            y = V*self._zs + (1-V)*y_bubble
            
            # Guess vapor flow rates
            self._v = y * V * self._mol

            # Solve
            self._vapor_mol[self._index] = self._solve_v(T, P)
            self._liquid_mol[self._index] = self._mol - self._v
            self.Q = stream.H/self._massnet
        
    def TV(self, T, V):
        stream = self._stream
        self.T = stream.T = T
        if self._N == 1: return self._TV_compound(T, V)
        P_dew, x_dew = self._dew_point.solve_Px(self._zs, T)
        P_bubble, y_bubble = self._bubble_point.solve_Py(self._zs, T)
        if V == 1:
            self._stream.P = P_dew
            self._vapor_mol[self._index] = self._mol
            self._liquid_mol[self._index] = 0
        elif V == 0:
            self._stream.P = P_bubble
            self._vapor_mol[self._index] = 0
            self._liquid_mol[self._index] = self._mol
        else:
            self.V = V
            self._v = (V*self._zs + (1-V)*y_bubble)*V*self._molnet
            self.P = stream.P = self.solver(self._V_at_P,
                                            P_bubble, P_dew, 0, 1,
                                            self.P, self.V,
                                            self.P_tol, self.V_tol)
            self._vapor_mol[self._index] = self._v
            self._liquid_mol[self._index] = self._mol - self._v
            self.Q = stream.H/self._massnet

    def TQ(self, T, Q):
        stream = self._stream
        stream.T = T
        if self._N == 1: return self._TQ_compound(T, Q)
        
        # Setup bounderies
        P_dew, x_dew = self._dew_point.solve_Px(self._zs, T)
        P_bubble, y_bubble = self._bubble_point.solve_Py(self._zs, T)
        index = self._index
        mol = self._mol
        vapor_mol = self._vapor_mol
        liquid_mol = self._liquid_mol
        Hin = Q + stream.H
        
        # Check if super heated vapor
        vapor_mol[index] = mol
        liquid_mol[index] = 0
        stream.P = P_dew
        H_dew = stream.H
        dH_dew = (Hin - H_dew)
        if dH_dew >= 0:
            stream.H = Hin
            return

        # Check if subcooled liquid
        vapor_mol[index] = 0
        liquid_mol[index] = mol
        stream.P = P_bubble
        H_bubble = stream.H
        dH_bubble = (Hin - H_bubble)
        if dH_bubble <= 0:
            stream.H = Hin
            return

        # Guess overall vapor fraction, and vapor flow rates
        V = self.V or dH_bubble/(H_dew - H_bubble)
        # Guess composition in the vapor is a weighted average of boiling points
        self._v = V*self._zs + (1-V)*y_bubble*V*self._molnet
        massnet = self._massnet
        self.Q = Hin/massnet
        self.P = stream.P = self.solver(self._Q_at_P,
                                        P_bubble, P_dew,
                                        H_bubble/massnet, H_dew/massnet,
                                        self.P, self.Q,
                                        self.P_tol, self.Q_tol) 
        self.T = T
    
    def PV(self, P, V):
        stream = self._stream
        stream.P = P
        if self._N == 1: return self._PV_compound(P, V)
        
        # Setup bounderies
        T_dew, x_dew = self._dew_point.solve_Tx(self._zs, P)
        T_bubble, y_bubble = self._bubble_point.solve_Ty(self._zs, P)
        
        index = self._index
        mol = self._mol
        vapor_mol = self._vapor_mol
        liquid_mol = self._liquid_mol
        
        if V == 1:
            stream.T = T_dew
            vapor_mol[index] = mol
            liquid_mol[index] = 0
            return
        elif V == 0:
            stream.T = T_bubble
            vapor_mol[index] = 0
            liquid_mol[index] = mol
            return
        else:
            self._v = (V*self._zs + (1-V)*y_bubble) * V * self._molnet
            self.V = V
            self.T = stream.T = self.solver(self._V_at_T,
                                            T_bubble, T_dew, 0, 1,
                                            self.T, V,
                                            self.T_tol, self.V_tol)
            vapor_mol[index] = self._v
            liquid_mol[index] = mol - self._v
            self.P = P
            self.Q = stream.H/self._massnet
    
    def PQ(self, P, Q):
        stream = self._stream
        stream.P = P
        if self._N == 1: return self._PQ_compound(P, Q)
        
        # Setup bounderies
        T_dew, x_dew = self._dew_point.solve_Tx(self._zs, P)
        T_bubble, y_bubble = self._bubble_point.solve_Ty(self._zs, P)
        
        index = self._index
        mol = self._mol
        vapor_mol = self._vapor_mol
        liquid_mol = self._liquid_mol
        Hin = Q + stream.H
        
        # Check if super heated vapor
        vapor_mol[index] = mol
        liquid_mol[index] = 0
        stream.T = T_dew
        H_dew = stream.H
        dH_dew = Hin - H_dew
        if dH_dew >= 0:
            stream.H = Hin
            return

        # Check if subcooled liquid
        vapor_mol[index] = 0
        liquid_mol[index] = mol
        stream.T = T_bubble
        H_bubble = stream.H
        dH_bubble = Hin - H_bubble
        if dH_bubble <= 0:
            stream.H = Hin
            return
        
        # Guess T, overall vapor fraction, and vapor flow rates
        self.V = V = self.V or dH_bubble/(H_dew - H_bubble)
        self._v = V*self._zs + (1-V)*y_bubble*V*self._molnet
        
        massnet = self._massnet
        self.Q = Hin/massnet
        self.T = stream.T = self.solver(self._Q_at_T,
                                        T_bubble, T_dew, 
                                        H_bubble/massnet, H_dew/massnet,
                                        self.T, self.Q,
                                        self.T_tol, self.Q_tol)
        self.P = P
    
    def _Q_at_T(self, T):
        self._stream.T = T
        self._vapor_mol[self._index] = self._solve_v(T, self._stream.P)
        self._liquid_mol[self._index] = self._mol - self._v
        return self._stream.H/self._massnet
    
    def _Q_at_P(self, P):
        self._stream.P = P
        self._vapor_mol[self._index] = self._solve_v(self._stream.T, P)
        self._liquid_mol[self._index] = self._mol - self._v
        return self._stream.H/self._massnet
    
    def _V_at_P(self, P):
        return self._solve_v(self._stream.T, P).sum()/self._molnet
    
    def _V_at_T(self, T):
        return self._solve_v(T, self._stream.P).sum()/self._molnet
    
    def _x_iter(self, x):
        self._Ks = self._Psat_P * self._gamma(x/x.sum(), self._stream.T)
        return self._zs/(1. + self._solve_V()*(self._Ks-1.))
    
    def _solve_v(self, T, P):
        """Solve for vapor mol"""
        gamma = self._gamma
        self._Psat_P = np.array([s.VaporPressure(T) for s in gamma.species])/P
        l = self._mol - self._v
        x = self.itersolver(self._x_iter, l/l.sum(), 1e-4)
        self._v = self._molnet*self.V*x/x.sum()*self._Ks            
        return self._v

    def _V_error(self, V):
        """Vapor fraction error."""
        return (self._zs*(self._Ks-1.)/(1.+V*(self._Ks-1.))).sum()

    def _solve_V_N(self):
        """Update V for N components."""
        self.V = self.solver(self._V_error, 0, 1,
                             self._V_error(0), self._V_error(1),
                             self.V, 0, 1e-4, 1e-7)
        return self.V

    def _solve_V_2(self):
        """Update V for 2 components."""
        z1, z2 = self._zs
        K1, K2 = self._Ks
        self.V = (-K1*z1 - K2*z2 + z1 + z2)/(K1*K2*z1 + K1*K2 *
                                             z2 - K1*z1 - K1*z2
                                             - K2*z1 - K2*z2 + z1 + z2)
        return self.V
    
    def _solve_V_3(self):
        """Update V for 3 components."""
        z1, z2, z3 = self._zs
        K1, K2, K3 = self._Ks
        self.V = ((-K1*K2*z1/2 - K1*K2*z2/2 - K1*K3*z1/2 - K1*K3*z3/2 + K1*z1
                   + K1*z2/2 + K1*z3/2 - K2*K3*z2/2 - K2*K3*z3/2 + K2*z1/2
                   + K2*z2 + K2*z3/2 + K3*z1/2 + K3*z2/2 + K3*z3 - z1 - z2 - z3
                   - (K1**2*K2**2*z1**2 + 2*K1**2*K2**2*z1*z2 + K1**2*K2**2*z2**2
                      - 2*K1**2*K2*K3*z1**2 - 2*K1**2*K2*K3*z1*z2 - 2*K1**2*K2*K3*z1*z3
                      + 2*K1**2*K2*K3*z2*z3 - 2*K1**2*K2*z1*z2 + 2*K1**2*K2*z1*z3
                      - 2*K1**2*K2*z2**2 - 2*K1**2*K2*z2*z3 + K1**2*K3**2*z1**2
                      + 2*K1**2*K3**2*z1*z3 + K1**2*K3**2*z3**2 + 2*K1**2*K3*z1*z2
                      - 2*K1**2*K3*z1*z3 - 2*K1**2*K3*z2*z3 - 2*K1**2*K3*z3**2
                      + K1**2*z2**2 + 2*K1**2*z2*z3 + K1**2*z3**2 - 2*K1*K2**2*K3*z1*z2
                      + 2*K1*K2**2*K3*z1*z3 - 2*K1*K2**2*K3*z2**2 - 2*K1*K2**2*K3*z2*z3
                      - 2*K1*K2**2*z1**2 - 2*K1*K2**2*z1*z2 - 2*K1*K2**2*z1*z3
                      + 2*K1*K2**2*z2*z3 + 2*K1*K2*K3**2*z1*z2 - 2*K1*K2*K3**2*z1*z3
                      - 2*K1*K2*K3**2*z2*z3 - 2*K1*K2*K3**2*z3**2 + 4*K1*K2*K3*z1**2
                      + 4*K1*K2*K3*z1*z2 + 4*K1*K2*K3*z1*z3 + 4*K1*K2*K3*z2**2
                      + 4*K1*K2*K3*z2*z3 + 4*K1*K2*K3*z3**2 + 2*K1*K2*z1*z2
                      - 2*K1*K2*z1*z3 - 2*K1*K2*z2*z3 - 2*K1*K2*z3**2 - 2*K1*K3**2*z1**2
                      - 2*K1*K3**2*z1*z2 - 2*K1*K3**2*z1*z3 + 2*K1*K3**2*z2*z3
                      - 2*K1*K3*z1*z2 + 2*K1*K3*z1*z3 - 2*K1*K3*z2**2 - 2*K1*K3*z2*z3
                      + K2**2*K3**2*z2**2 + 2*K2**2*K3**2*z2*z3 + K2**2*K3**2*z3**2
                      + 2*K2**2*K3*z1*z2 - 2*K2**2*K3*z1*z3 - 2*K2**2*K3*z2*z3
                      - 2*K2**2*K3*z3**2 + K2**2*z1**2 + 2*K2**2*z1*z3
                      + K2**2*z3**2 - 2*K2*K3**2*z1*z2 + 2*K2*K3**2*z1*z3
                      - 2*K2*K3**2*z2**2 - 2*K2*K3**2*z2*z3 - 2*K2*K3*z1**2
                      - 2*K2*K3*z1*z2 - 2*K2*K3*z1*z3 + 2*K2*K3*z2*z3 + K3**2*z1**2
                      + 2*K3**2*z1*z2 + K3**2*z2**2)**0.5/2)
                   / (K1*K2*K3*z1 + K1*K2*K3*z2 + K1*K2*K3*z3 - K1*K2*z1
                      - K1*K2*z2 - K1*K2*z3 - K1*K3*z1 - K1*K3*z2 - K1*K3*z3
                      + K1*z1 + K1*z2 + K1*z3 - K2*K3*z1 - K2*K3*z2 - K2*K3*z3
                      + K2*z1 + K2*z2 + K2*z3 + K3*z1 + K3*z2 + K3*z3 - z1 - z2 - z3))
        return self.V

    def __repr__(self):
        return f"{repr(self._stream)}.VLE"
