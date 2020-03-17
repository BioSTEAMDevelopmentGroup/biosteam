# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 02:11:08 2020

@author: yoelr
"""
__all__ = ('ExponentialFunctor',)

class ExponentialFunctor:
    """
    Create an ExponentialFunctor object for computing equations of the 
    form :math:`f(S) = A \cdot S^n`.
    
    Parameters
    ----------
    A : float
        Linear coefficient.
    n : float
        Exponential coefficient.
    
    Attributes
    ----------
    A : float
        Linear coefficient.
    n : float
        Exponential coefficient.
        
    Examples
    --------
    >>> from biosteam.utils import ExponentialFunctor
    >>> f_exp = ExponentialFunctor(A=10.0, n=2.0)
    >>> f_exp(2.0)
    40.0
    """    
    __slots__ = ('A', 'n')

    def __init__(self, A, n):
        self.A = A
        self.n = n

    def __call__(self, S):
        return self.A * S ** self.n

    def __repr__(self):
        return f"{type(self).__name__}(A={self.A}, n={self.n})"
    
