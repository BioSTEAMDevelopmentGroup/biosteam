#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 08:58:41 2018

@author: Yoel Rene Cortes-Pena
"""
from math import log, exp
from ..thermo.unifac import DOUFSG, DOUFIP2016

__all__ = ('Dortmund',)


# %% Activity Coefficients

class Dortmund:
    __slots__ = ('_species', 'chemgroups', 'rs', 'qs', 'groupcounts')
    subgroups = DOUFSG
    interactions = DOUFIP2016
    cached = {}
    
    def __init__(self, *species):
        self._species = self.chemgroups = self.rs = self.qs = self.groupcounts = None
        if species: self.species = species
    
    @property
    def species(self):
        return self._species
    
    @species.setter
    def species(self, species):
        if self._species != species:
            self._species = species
            species = tuple(species)
            if species in self.cached:
                self.chemgroups, self.rs, self.qs, self.groupcounts = self.cached[species]
            else:
                self._load_species(species)
                
    
    def __call__(self, xs, T):
        """Return UNIFAC coefficients.
        
        Parameters
        ----------
        xs : array_like
            Molar fractions
        T : float
            Temperature (K)
        
        """
        groupcounts = self.groupcounts
        gckeys = tuple(groupcounts)
        xs_chemgroups = tuple(zip(xs, self.chemgroups))
        # Sum the denominator for calculating Xs
        sum_ = sum
        group_sum = sum_([count*x for x, g in xs_chemgroups for count in g.values()])
    
        # Caclulate each numerator for calculating Xs
        group_count_xs = {}
        for group in gckeys:
            tot_numerator = sum_([x*g[group] for x, g in xs_chemgroups if group in g])
            group_count_xs[group] = tot_numerator/group_sum
    
        loggammacs = self.loggammacs(self.qs, self.rs, xs)
        subgroups = self.subgroups
        Q_sum_term = sum_([subgroups[group].Q*group_count_xs[group] for group in gckeys])
        area_fractions = {group: subgroups[group].Q*group_count_xs[group]/Q_sum_term
                          for group in gckeys}
    
        psi = self.psi
        interactions = self.interactions
        UNIFAC_psis = {k: {m: (psi(T, m, k, subgroups, interactions))
                           for m in gckeys} for k in gckeys}
    
        loggamma_groups = {}
        for k in gckeys:
            sum1, sum2 = 0., 0.
            for m in gckeys:
                sum1 += area_fractions[m]*UNIFAC_psis[k][m]
                sum3 = sum_([area_fractions[n]*UNIFAC_psis[m][n]
                            for n in groupcounts])
                sum2 -= area_fractions[m]*UNIFAC_psis[m][k]/sum3
            loggamma_groups[k] = subgroups[k].Q*(1. - log(sum1) + sum2)
    
        loggammars = []
        for groups in self.chemgroups:
            gkeys = groups.keys()
            chem_group_sum = sum_(groups.values())
            chem_group_count_xs = {group: count/chem_group_sum
                                   for group, count in groups.items()}
    
            Q_sum_term = sum_([subgroups[group].Q*chem_group_count_xs[group]
                               for group in gkeys])
            chem_area_fractions = {group: subgroups[group].Q*chem_group_count_xs[group]/Q_sum_term
                                   for group in gkeys}
            chem_loggamma_groups = {}
            for k in gkeys:
                sum1, sum2 = 0., 0.
                for m in groups:
                    sum1 += chem_area_fractions[m]*UNIFAC_psis[k][m]
                    sum3 = sum_([chem_area_fractions[n]
                                * UNIFAC_psis[m][n] for n in groups])
                    sum2 -= chem_area_fractions[m]*UNIFAC_psis[m][k]/sum3
    
                chem_loggamma_groups[k] = subgroups[k].Q*(1. - log(sum1) + sum2)
    
            tot = sum_([count*(loggamma_groups[group] - chem_loggamma_groups[group])
                       for group, count in groups.items()])
            loggammars.append(tot)
    
        return [exp(sum_(ij)) for ij in zip(loggammacs, loggammars)]
    
    def _load_species(self, species):
        self.chemgroups = chemgroups = [s.UNIFAC_Dortmund_groups for s in species]
        self.rs = rs = []
        self.qs = qs = []
        self.groupcounts = groupcounts = {}
        subgroups = self.subgroups
        for groups in chemgroups:
            ri = 0.
            qi = 0.
            for group, count in groups.items():
                ri += subgroups[group].R*count
                qi += subgroups[group].Q*count
                if group in groupcounts: groupcounts[group] += count
                else: groupcounts[group] = count
            rs.append(ri)
            qs.append(qi)
        self.cached[species] = (chemgroups, rs, qs, groupcounts)

    def loggammacs(self, qs, rs, xs):
        rsxs = 0.; qsxs = 0.; rsxs2 = 0.; rs_34 = []; loggammacs = []
        for xi, ri, qi in zip(xs, rs, qs):
            rsxs += ri*xi
            qsxs += qi*xi
            ri_34 = ri**0.75
            rs_34.append(ri_34)
            rsxs2 += ri_34*xi
        for qi, ri, ri_34 in zip(qs, rs, rs_34):
            Vi = ri/rsxs
            Fi = qi/qsxs
            Vi2 = ri_34/rsxs2
            loggammacs.append(1. - Vi2 + log(Vi2) - 5.*qi*(1. - Vi/Fi + log(Vi/Fi)))
        return loggammacs
    
    def psi(self, T, subgroup1, subgroup2, subgroup_data, interaction_data):
        try: a, b, c = interaction_data[subgroup_data[subgroup1].main_group_id] \
                                       [subgroup_data[subgroup2].main_group_id]
        except: return 1.
        return exp((-a/T - b - c*T))
    
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



