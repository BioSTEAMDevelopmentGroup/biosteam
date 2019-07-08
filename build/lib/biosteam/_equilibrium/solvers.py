# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 18:40:05 2019

@author: yoelr
"""
from .._utils import isbetween

__all__ = ('VLEsolver',)

# def boundsolve(f, x0, x, x1, y0, y, y1, yval, xtol, ytol):
#     if y1 < 0: x0, y0, x1, y1 = x1, y1, x0, y0
#     if not (isbetween(x0, x, x1) and (yval-10*ytol < y < yval+10*ytol)):
#         x = x0 + (yval-y0)*(x1-x0)/(y1-y0)
#     y = f(x)
#     x_ = x0
#     while abs(x-x_) > xtol and abs(yval-y) > ytol:
#         if y > yval:
#             x_ = x1 = x
#             y1 = y
#         else:
#             x_ = x0 = x
#             y0 = y
#         dx = x1-x0
#         x = x0 + (yval-y0)*(dx)/(y1-y0)
#         y = f(x)
#     return x, y

class VLEsolver:
    __slots__ = ('T', 'P', 'Q', 'V')
    
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
        while abs(x-x_) > xtol and abs(yval-y) > ytol:
            if y > yval:
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
