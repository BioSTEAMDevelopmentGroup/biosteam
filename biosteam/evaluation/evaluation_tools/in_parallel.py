# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import multiprocessing as mp
import pandas as pd

__all__ = ('evaluate_coordinate_in_parallel',)

def evaluate_coordinate_in_parallel(f_evaluate_at_coordinate, coordinate,
                                    metrics, multi_coordinate=False,
                                    name=None, names=None, xlfile=None): # pragma: no cover
    # In parallel
    pool = mp.Pool(mp.cpu_count())
    map = pool.starmap if multi_coordinate else pool.map
    tables = map(f_evaluate_at_coordinate, coordinate)
    pool.close()
    
    # Initialize data containers
    table = tables[0]
    N_samples, N_parameters = table.shape
    N_coordinate = len(coordinate)
    metric_data = {i.index: np.zeros([N_samples, N_coordinate])
                   for i in metrics}
    try:
        for n, table in enumerate(tables):
            for metric in metric_data:
                metric_data[metric][:, n] = table[metric]
    except:
        return tables, metric_data
    
    if xlfile:
        if multi_coordinate:
            columns = pd.MultiIndex.from_tuples(coordinate,
                                                names=names or name)
        else:
            columns = pd.Index(coordinate, name=name or names)
        
        # Save data to excel
        data = pd.DataFrame(data=np.zeros([N_samples, N_coordinate]),
                            columns=columns)
        
        with pd.ExcelWriter(xlfile) as writer:
            for i, metric in zip(metrics, metric_data):
                data[:] = metric_data[metric]
                data.to_excel(writer, sheet_name=i.name)
    
    return metric_data
        