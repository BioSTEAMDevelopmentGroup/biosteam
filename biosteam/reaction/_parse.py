# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 20:47:24 2019

@author: yoelr
"""

__all__ = ('str2dct', 'dct2str', 'dct2arr', 'arr2dct',  'str2arr', 'arr2str')

def str2dct(reaction) -> dict:
    left, right = reaction.split('->')
    reactants = left.split('+')
    products = right.split('+')
    dct = {}
    for nID in reactants:
        for i, letter in enumerate(nID):
            if letter.isalpha(): break
        if i: dct[nID[i:]] = -float(nID[:i])
        else: dct[nID] = -1
    for nID in products:
        for i, letter in enumerate(nID):
            if letter.isalpha(): break
        if i: dct[nID[i:]] = float(nID[:i])
        else: dct[nID] = 1
    return dct

def dct2str(dct):
    left = []
    right = []
    for ID, N in dct.items():
        N_int = int(N)
        if N_int == N: N = N_int
        if N == -1: left.append(ID)
        elif N == 1: right.append(ID)
        elif N < 0: left.append(f"{-N} {ID}")
        else: right.append(f"{N} {ID}")
    left = ' + '.join(left)
    right = ' + '.join(right)
    reaction = left + ' -> ' + right
    return reaction

def dct2arr(dct, species):
    return species.array(dct.keys(), [*dct.values()])

def arr2dct(arr, species):
    dct = {}
    for N, ID in zip(arr, species._IDs):
        if N: dct[ID] = N
    return dct

def str2arr(reaction, species):
    return dct2arr(str2dct(reaction), species)

def arr2str(arr, species):
    return dct2str(arr2dct(arr, species))
