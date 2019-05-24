# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:13:50 2019

Functions for inspecting separations, including functions to calculate component splits, and partition coefficients between streams

@author: yoelr
"""

import pandas as pd
import numpy as np

__all__ = ('compound_split', 'partition_coefficient')

def compound_split(streams):
    """Return compound split fractions for each stream in outs as a pandas DataFrame."""
    mol_out = sum(s._mol for s in streams)
    mol_out[mol_out==0] = 1
    return pd.DataFrame(data=np.array([i.mol/mol_out for i in streams]).transpose(),
                        index=streams[0]._species._IDs,
                        columns=[i.ID for i in streams])

def partition_coefficient(streams):
    """Return partition coefficient information as a pandas Series."""
    ph1, ph2 = streams

    # Account for light and heave keys
    ph1zero = ph1.mol == 0
    ph2zero = ph2.mol == 0
    HNK = (ph1zero) & (~ph2zero)
    LNK = (ph2zero) & (~ph1zero)

    # Do not include heavy or light non keys
    pos = ~ph1zero & ~ph2zero
    ph1_mol = ph1.mol[pos]
    ph2_mol = ph2.mol[pos]
    ph1_molnet = ph1_mol.sum()
    ph2_molnet = ph2_mol.sum()
    ph1_molfrac = ph1_mol/ph1_molnet
    ph2_molfrac = ph2_mol/ph2_molnet
    IDs = np.array(ph1._species._IDs)
    data = {} 
    data['Top fraction'] = ph1_molnet/(ph1_molnet + ph2_molnet)
    data.update(zip(IDs[pos], ph1_molfrac/ph2_molfrac))
    data['Light non-keys'] = IDs[LNK]
    data['Heavy non-keys'] = IDs[HNK]
    return pd.Series(data)