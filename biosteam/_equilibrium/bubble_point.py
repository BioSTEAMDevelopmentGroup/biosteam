# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 21:30:33 2019

@author: yoelr
"""
from numpy import asarray, array
from flexsolve import wegstein_secant, aitken_secant

__all__ = ('BubblePoint',)

class BubblePoint:
    __slots__ = ('gamma', 'P', 'T', 'y')
    rootsolver = staticmethod(aitken_secant)
    
    def __init__(self, gamma):
        self.gamma = gamma
    
    def _T_error(self, T, z_over_P, z, VPs):
        self.y =  z_over_P * array([i(T) for i in VPs]) * self.gamma(z, T)
        return 1. - self.y.sum()
    
    def _P_error(self, P, T, z, Psat_gamma):
        self.y = z * Psat_gamma / P
        return 1. - self.y.sum()
    
    def solve_Ty(self, z, P):
        """Bubble point at given composition and pressure

        Parameters
        ----------
        z : array_like
            Liquid phase composition.
        P : float
            Pressure (Pa).
        
        Returns
        -------
        T : float 
            Bubble point temperature (K)
        y : numpy.ndarray
            Composition of the vapor phase.

        Examples
        --------
        >>> from biosteam import Species, BubblePoint, Dortmund
        >>> gamma = Dortmund(*Species('Ethanol', 'Water'))
        >>> bp = BubblePoint(gamma)
        >>> bp.solve_Ty(z=(0.6, 0.4), P=101325)
        (352.2820850833474, array([0.703, 0.297]))
        
        """
        z = asarray(z)
        self.P = P
        try:
            self.T = self.rootsolver(self._T_error, self.T, self.T+0.01, 1e-6,
                                     args=(z/P, z, [s.VaporPressure for s in self.gamma.species]))
        except:
            T = (z * [s.Tb for s in self.gamma.species]).sum()
            self.T = self.rootsolver(self._T_error, T, T+0.01, 1e-6,
                                     args=(z/P, z, [s.VaporPressure for s in self.gamma.species]))
        self.y = self.y/self.y.sum()
        return self.T, self.y
    
    def solve_Py(self, z, T):
        """Bubble point at given composition and temperature.

        Parameters
        ----------
        z : array_like
            Liquid phase composotion.
        T : float
            Temperature (K).
        
        Returns
        -------
        P : float
            Bubble point pressure (Pa).
        y : numpy.ndarray
            Vapor phase composition.

        Examples
        --------
        >>> from biosteam import Species, BubblePoint, Dortmund
        >>> gamma = Dortmund(*Species('Ethanol', 'Water'))
        >>> bp = BubblePoint(gamma)
        >>> bp.solve_Py(z=(0.703, 0.297), T=352.28)
        (103494.17209657285, array([0.757, 0.243]))
        
        """
        z = asarray(z)
        Psat = array([i.VaporPressure(T) for i in self.gamma.species])
        Psat_gamma =  Psat * self.gamma(z, T)
        self.T = T
        try:
            self.P = self.rootsolver(self._P_error, self.P, self.P-1, 1e-2,
                                     args=(T, z, Psat_gamma))
        except:
            P = (z * Psat).sum()
            self.P = self.rootsolver(self._P_error, P, P-1, 1e-2,
                                     args=(T, z, Psat_gamma))
        self.y = self.y/self.y.sum()
        return self.P, self.y
    
    def __repr__(self):
        return f"<{type(self).__name__}: gamma={self.gamma}>"
    

    