# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:40:05 2019

@author: yoelr
"""
from .._utils import isbetween, boundsolve, wegstein
from .._exceptions import SolverError
import numpy as np

__all__ = ('VLEsolver', 'solve_v', 'V_2N', 'V_3N', 'V_error')

def V_2N(zs, Ks):
    """Solution for 2 component flash vessel."""
    z1, z2 = zs
    K1, K2 = Ks
    return (-K1*z1 - K2*z2 + z1 + z2)/(K1*K2*z1 + K1*K2 *
                                       z2 - K1*z1 - K1*z2
                                       - K2*z1 - K2*z2 + z1 + z2)
    
def V_3N(zs, Ks):
    """Solution for 3 component flash vessel."""
    z1, z2, z3 = zs
    K1, K2, K3 = Ks
    return (-K1*K2*z1/2 - K1*K2*z2/2 - K1*K3*z1/2 - K1*K3*z3/2 + K1*z1 + K1*z2/2 + K1*z3/2 - K2*K3*z2/2 - K2*K3*z3/2 + K2*z1/2 + K2*z2 + K2*z3/2 + K3*z1/2 + K3*z2/2 + K3*z3 - z1 - z2 - z3 - (K1**2*K2**2*z1**2 + 2*K1**2*K2**2*z1*z2 + K1**2*K2**2*z2**2 - 2*K1**2*K2*K3*z1**2 - 2*K1**2*K2*K3*z1*z2 - 2*K1**2*K2*K3*z1*z3 + 2*K1**2*K2*K3*z2*z3 - 2*K1**2*K2*z1*z2 + 2*K1**2*K2*z1*z3 - 2*K1**2*K2*z2**2 - 2*K1**2*K2*z2*z3 + K1**2*K3**2*z1**2 + 2*K1**2*K3**2*z1*z3 + K1**2*K3**2*z3**2 + 2*K1**2*K3*z1*z2 - 2*K1**2*K3*z1*z3 - 2*K1**2*K3*z2*z3 - 2*K1**2*K3*z3**2 + K1**2*z2**2 + 2*K1**2*z2*z3 + K1**2*z3**2 - 2*K1*K2**2*K3*z1*z2 + 2*K1*K2**2*K3*z1*z3 - 2*K1*K2**2*K3*z2**2 - 2*K1*K2**2*K3*z2*z3 - 2*K1*K2**2*z1**2 - 2*K1*K2**2*z1*z2 - 2*K1*K2**2*z1*z3 + 2*K1*K2**2*z2*z3 + 2*K1*K2*K3**2*z1*z2 - 2*K1*K2*K3**2*z1*z3 - 2*K1*K2*K3**2*z2*z3 - 2*K1*K2*K3**2*z3**2 + 4*K1*K2*K3*z1**2 + 4*K1*K2*K3*z1*z2 + 4*K1*K2*K3*z1*z3 + 4*K1*K2*K3*z2**2 + 4*K1*K2*K3*z2*z3 + 4*K1*K2*K3*z3**2 + 2*K1*K2*z1*z2 - 2*K1*K2*z1*z3 - 2*K1*K2*z2*z3 - 2*K1*K2*z3**2 - 2*K1*K3**2*z1**2 - 2*K1*K3**2*z1*z2 - 2*K1*K3**2*z1*z3 + 2*K1*K3**2*z2*z3 - 2*K1*K3*z1*z2 + 2*K1*K3*z1*z3 - 2*K1*K3*z2**2 - 2*K1*K3*z2*z3 + K2**2*K3**2*z2**2 + 2*K2**2*K3**2*z2*z3 + K2**2*K3**2*z3**2 + 2*K2**2*K3*z1*z2 - 2*K2**2*K3*z1*z3 - 2*K2**2*K3*z2*z3 - 2*K2**2*K3*z3**2 + K2**2*z1**2 + 2*K2**2*z1*z3 + K2**2*z3**2 - 2*K2*K3**2*z1*z2 + 2*K2*K3**2*z1*z3 - 2*K2*K3**2*z2**2 - 2*K2*K3**2*z2*z3 - 2*K2*K3*z1**2 - 2*K2*K3*z1*z2 - 2*K2*K3*z1*z3 + 2*K2*K3*z2*z3 + K3**2*z1**2 + 2*K3**2*z1*z2 + K3**2*z2**2)**0.5/2)/(K1*K2*K3*z1 + K1*K2*K3*z2 + K1*K2*K3*z3 - K1*K2*z1 - K1*K2*z2 - K1*K2*z3 - K1*K3*z1 - K1*K3*z2 - K1*K3*z3 + K1*z1 + K1*z2 + K1*z3 - K2*K3*z1 - K2*K3*z2 - K2*K3*z3 + K2*z1 + K2*z2 + K2*z3 + K3*z1 + K3*z2 + K3*z3 - z1 - z2 - z3)

def V_error(V, zs, Ks):
    """Vapor fraction error"""
    return (zs*(Ks-1.)/(1.+V*(Ks-1.))).sum()

def solve_v(v, T, P, mol, molnet, zs, N, species, f_activity):
    """Solve for vapor mol"""
    Psat_P = np.array([s.VaporPressure(T) for s in species])/P
    V = v.sum()/molnet
    l = mol - v
    x = l/l.sum()
    if N == 2:
        solve_V = V_2N
    elif N == 3:
        solve_V = V_3N
    else:
        solve_V = lambda zs, Ks: boundsolve(V_error, 0, V, 1, 0.0001, 0.001, args=(zs, Ks))
    Ks = None
    def f(x):
        nonlocal Ks, V
        Ks = Psat_P * f_activity(species, x/x.sum(), T)
        V = solve_V(zs, Ks)
        return zs/(1. + V*(Ks-1.))
    x = wegstein(f, x, 0.001)
    return molnet*V*x/x.sum()*Ks

class VLEsolver:
    """Create a VLEsolver object for solving VLE."""
    __slots__ = ('T', 'P', 'Q', 'V')
    
    maxiter = 30
    tolerance = {'T': 0.005,
                 'P': 10,
                 'Q': 1,
                 'V': 0.0005}
    
    def __init__(self):
        self.T = self.P = self.Q = self.V = 0
    
    def __call__(self, xvar, yvar, f, x0, x1, y0, y1, yval):
        xtol = self.tolerance[xvar]
        ytol = self.tolerance[yvar]
        x = getattr(self, xvar)
        y = getattr(self, yvar)
        if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
        if not (isbetween(x0, x, x1) and (yval-10*ytol < y < yval+10*ytol)):
            x = x0 + (yval-y0)*(x1-x0)/(y1-y0)
        y = f(x)
        x_ = x0
        it = 0
        maxiter = self.maxiter
        while abs(x-x_) > xtol and abs(yval-y) > ytol:
            it += 1
            if it > maxiter:
                raise SolverError('failed to converge')
            elif y > yval:
                x_ = x1 = x
                y1 = y
            else:
                x_ = x0 = x
                y0 = y
            dx = x1-x0
            x = x0 + (yval-y0)*(dx)/(y1-y0)
            y = f(x)
        setattr(self, xvar, x)
        setattr(self, yvar, y)
        return x

    def __repr__(self):
        return f"<{type(self).__name__}: T={self.T:.2f} K, P={int(self.P)} Pa, Q={self.Q:.3g} kJ/kg, V={self.V:.2f}>"
