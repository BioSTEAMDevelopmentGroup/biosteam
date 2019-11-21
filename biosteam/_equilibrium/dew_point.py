# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 22:15:30 2019

@author: yoelr
"""

from numpy import asarray, array
from flexsolve import wegstein_secant, aitken_secant, wegstein, aitken

__all__ = ('DewPoint',)

class DewPoint:
    __slots__ = ('gamma', 'P', 'T', 'x')
    rootsolver = staticmethod(aitken_secant)
    itersolver = staticmethod(aitken)
    def __init__(self, gamma):
        self.gamma = gamma
    
    def _x_iter_at_P(self, x, T, zP, VPs):
        return zP/(array([i(T) for i in VPs]) * self.gamma(x/x.sum(), T))
    
    def _x_iter_at_T(self, x, T, P, z, Psat):
        return z*P/(Psat * self.gamma(x/x.sum(), T))
    
    def _T_error(self, T, zP, VPs):
        self.x = self.itersolver(self._x_iter_at_P, self.x, 1e-5, args=(T, zP, VPs))
        return 1 - self.x.sum()
    
    def _P_error(self, P, T, z, Psat):
        self.x = self.itersolver(self._x_iter_at_T, self.x, 1e-5, args=(T, P, z, Psat))
        return 1 - self.x.sum()
    
    def solve_Tx(self, z, P):
        """Dew point given composition and pressure.

        Parameters
        ----------
        y : array_like
            Vapor phase composition.

        P : float
            Pressure (Pa).

        Returns
        -------
        T : float
            Dew point temperature (K).
        x : numpy.ndarray
            Liquid phase composition.

        Examples
        --------
        >>> from biosteam import Species, DewPoint, Dortmund
        >>> gamma = Dortmund(*Species('Ethanol', 'Water'))
        >>> dp = DewPoint(gamma)
        >>> dp.solve_Tx(z=(0.5, 0.5), P=101325)
        (357.45184742263075, array([0.151, 0.849]))
        """
        z = asarray(z)
        self.P = P
        try:
            self.T = self.rootsolver(self._T_error, self.T, self.T-0.01, 1e-6,
                                     args=(P*z, [s.VaporPressure for s in self.gamma.species]))
        except:
            self.x = z.copy()
            T = (z * [s.Tb for s in self.gamma.species]).sum()
            try:
                self.T = self.rootsolver(self._T_error, T, T-0.01, 1e-6,
                                         args=(P*z, [s.VaporPressure for s in self.gamma.species]))
            except:
                self.x = z.copy()
                T_guess = min([s.Tb for s in self.gamma.species])
                try:
                    self.T = self.rootsolver(self._T_error, T_guess, T_guess-0.01, 1e-6,
                                             args=(P*z, [s.VaporPressure for s in self.gamma.species]))
                except:
                    self.T = T
                    self.x = z.copy()
                
        self.x = self.x/self.x.sum()
        return self.T, self.x
    
    def solve_Px(self, z, T):
        """Dew point given composition and temperature.

        Parameters
        ----------
        y : array_like
            Vapor phase composition.
        T : float
            Temperature (K).
        
        Returns
        -------
        P : float
            Dew point pressure (Pa).
        x : numpy.ndarray
            Liquid phase composition.

        Examples
        --------
        >>> from biosteam import Species, DewPoint, Dortmund
        >>> gamma = Dortmund(*Species('Ethanol', 'Water'))
        >>> dp = DewPoint(gamma)
        >>> dp.solve_Px(z=(0.703, 0.297), T=352.28)
        (111366.15384513882, array([0.6, 0.4]))
 
       """
        z = asarray(z)
        Psat = array([i.VaporPressure(T) for i in self.gamma.species])
        self.T = T
        try:
            self.P = self.rootsolver(self._P_error, self.P, self.P+1, 1e-2,
                                     args=(T, z, Psat))
        except:
            P = (z * Psat).sum()
            self.x = z.copy()
            self.P = self.rootsolver(self._P_error, P, P+1, 1e-2,
                                     args=(T, z, Psat))
        self.x = self.x/self.x.sum()
        return self.P, self.x
    
    def __repr__(self):
        return f"<{type(self).__name__}: gamma={self.gamma}>"