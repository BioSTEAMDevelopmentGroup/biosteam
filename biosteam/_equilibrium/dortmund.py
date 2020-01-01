#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 08:58:41 2018

@author: Yoel Rene Cortes-Pena
"""
import numpy as np
from numba import njit
from ..thermo.unifac import DOUFSG, DOUFIP2016

__all__ = ('Dortmund',)

# %% Utilities

def chemgroup_array(chemgroups, index):
    M = len(chemgroups)
    N = len(index)
    array = np.zeros((M, N))
    for i, groups in enumerate(chemgroups):
        for group, count in groups.items():
            array[i, index[group]] = count
    return array
@njit
def group_activity_coefficients(x, chemgroups, loggammacs,
                                Qs, psis, cQfs, gpsis):
    weighted_counts = chemgroups.transpose() @ x
    Q_fractions = Qs * weighted_counts 
    Q_fractions /= Q_fractions.sum()
    Q_psis = psis * Q_fractions
    sum1 = Q_psis.sum(1)
    sum2 = -(psis.transpose() / sum1) @ Q_fractions
    loggamma_groups = Qs * (1. - np.log(sum1) + sum2)
    sum1 = cQfs @ gpsis.transpose()
    sum1 = np.where(sum1==0, 1., sum1)
    fracs = - cQfs / sum1
    sum2 = fracs @ gpsis
    chem_loggamma_groups = Qs*(1. - np.log(sum1) + sum2)
    loggammars = ((loggamma_groups - chem_loggamma_groups) * chemgroups).sum(1)
    return np.exp(loggammacs + loggammars)

def get_interaction(all_interactions, i, j, no_interaction):
    if i==j:
        return no_interaction
    try:
        return all_interactions[i][j]
    except:
        return no_interaction

# %% Activity Coefficients

class Dortmund:
    __slots__ = ('_rs', '_qs', '_Qs','_chemgroups',
                 '_group_psis',  '_chem_Qfractions',
                 '_group_mask', '_interactions',
                 '_species')
    all_subgroups = DOUFSG
    all_interactions = DOUFIP2016
    _no_interaction = (0., 0., 0.)
    cache = {}
    
    def __init__(self, *species):
        self._species = None
        if species: self.species = species
        
    @property
    def species(self):
        return self._species
    
    @species.setter
    def species(self, species):
        if self._species != species:
            self._species = species
            species = tuple(species)
            if species in self.cache:
                other = self.cache[species]
                for i in self.__slots__:
                    setattr(self, i, getattr(other, i))
            else:
                self._load_species(species)
           
    def __call__(self, x, T):
        """Return UNIFAC coefficients.
        
        Parameters
        ----------
        x : array_like
            Molar fractions
        T : float
            Temperature (K)
        
        """
        x = np.asarray(x)
        psis = self.psi(T, self._interactions.copy())
        self._group_psis[self._group_mask] =  psis[self._group_mask]
        return group_activity_coefficients(x, self._chemgroups,
                                           self.loggammacs(self._qs, self._rs, x),
                                           self._Qs, psis,
                                           self._chem_Qfractions,
                                           self._group_psis)
    
    @staticmethod
    @njit
    def loggammacs(qs, rs, x):
        r_net = (x*rs).sum()
        q_net = (x*qs).sum()
        rs_p = rs**0.75
        r_pnet = (rs_p*x).sum()
        Vs = rs/r_net
        Fs = qs/q_net
        Vs_over_Fs = Vs/Fs
        Vs_p = rs_p/r_pnet
        return 1. - Vs_p + np.log(Vs_p) - 5.*qs*(1. - Vs_over_Fs + np.log(Vs_over_Fs))
    
    @staticmethod
    @njit
    def psi(T, abc):
        abc[:, :, 0] /= T
        abc[:, :, 2] *= T
        return np.exp(-abc.sum(2)) 
    
    def _load_species(self, species):
        chemgroups = [s.UNIFAC_Dortmund_groups for s in species]
        all_groups = set()
        for groups in chemgroups: all_groups.update(groups)
        index = {group:i for i,group in enumerate(all_groups)}
        chemgroups = chemgroup_array(chemgroups, index)
        all_subgroups = self.all_subgroups
        subgroups = [all_subgroups[i] for i in all_groups]
        main_group_ids = [i.main_group_id for i in subgroups]
        self._Qs = Qs = np.array([i.Q for i in subgroups])
        Rs = np.array([i.R for i in subgroups])
        self._rs = chemgroups @ Rs
        self._qs = chemgroups @ Qs
        self._chemgroups = chemgroups
        chem_Qs = Qs * chemgroups
        self._chem_Qfractions = cQfs = chem_Qs/chem_Qs.sum(1, keepdims=True)
        all_interactions = self.all_interactions
        N_groups = len(all_groups)
        group_shape = (N_groups, N_groups)
        self._interactions = np.array([[get_interaction(all_interactions, i, j, self._no_interaction)
                                        for i in main_group_ids]
                                       for j in main_group_ids])
        # Psis array with only symmetrically available groups
        self._group_psis = np.zeros(group_shape, dtype=float)
        # Make mask for retrieving symmetrically available groups
        rowindex = np.arange(N_groups, dtype=int)
        indices = [rowindex[rowmask] for rowmask in cQfs != 0]
        self._group_mask = group_mask = np.zeros(group_shape, dtype=bool)
        for index in indices:
            for i in index:
                group_mask[i, index] = True
        self._species = species

    
    def __repr__(self):
        return f"{type(self).__name__}({', '.join([i.ID for i in self.species])})"
    
    
# def UNIFAC(self, xs, T):
#     return UNIFAC_Coeffictients(self, xs, T, UFSG, UFIP, UNIFAC_psi, loggammacs_UNIFAC)

# def UNIFAC_LL(self, xs, T):
#     """For LLE"""
#     return UNIFAC_Coeffictients(self, xs, T, UFSG, UFLLIP, UNIFAC_psi, loggammacs_UNIFAC)

# def loggammacs_UNIFAC(qs, rs, xs):
#     rsxs = sum([ri*xi for ri, xi in zip(rs, xs)])
#     Vis = [ri/rsxs for ri in rs]
#     qsxs = sum([qi*xi for qi, xi in zip(qs, xs)])
#     Fis = [qi/qsxs for qi in qs]

#     loggammacs = [1. - Visi + log(Visi) - 5.*qsi*(1. - Visi/Fisi + log(Visi/Fisi))
#                   for Visi, Fisi, qsi in zip(Vis, Fis, qs)]
#     return loggammacs

# def UNIFAC_psi(T, subgroup1, subgroup2, subgroup_data, interaction_data):
#     try:
#         return exp(-interaction_data[subgroup_data[subgroup1].main_group_id] \
#                                     [subgroup_data[subgroup2].main_group_id]/T)
#     except:
#         return 1



