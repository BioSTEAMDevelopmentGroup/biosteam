# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import matplotlib.pyplot as plt
from .table import cost_table
import numpy as np

__all__ = ['plot_cost_summary']

def plot_cost_summary(system, cost='Capital'): # pragma: no coverage
    cost = cost.capitalize()
    if cost == 'Capital':
        colkey = 'Fixed Capital Investment (10^6 USD)'
        ylabel = 'Fixed Capital Investment\n[$10^6$ USD]'
    elif cost == 'Utility':
        colkey = 'Utility Cost (10^6 USD/yr)'
        ylabel = 'Utility Cost\n[$10^6$ USD]'
    else:
        raise ValueError(f"Argument, cost, must be either 'Capital' or 'Utility'")
    # Get data
    df = cost_table(system)
    Type = df['Type'] # Column of unit types
    C_unit = df[colkey] # Column of unit costs
    C_Type = [] # Cost of all units of a type
    t_old = Type[0]
    r = 0
    c_type = 0
    Type_ = [] # All types
    for t in Type:
        if t == t_old:
            c_type += C_unit[r]
        else:
            C_Type.append(c_type)
            Type_.append(t_old)
            c_type = C_unit[r]
            t_old = t
        r += 1
    if t != Type_[-1]:
        C_Type.append(c_type)
        Type_.append(t)
        
    # Make graph
    xpos = np.arange(len(Type_))*1.2
    max_ = max(C_Type)
    ypos = np.arange(max_*1.2)
    plt.figure(figsize=(18,5))
    plt.bar(xpos, C_Type, align='center', alpha=0.5)
    plt.xticks(xpos, Type_, fontsize='11', rotation=30, ha='right')
    plt.yticks(ypos, fontsize='11')
    plt.xlabel('Unit Operation Type', fontsize='14')
    plt.ylabel(ylabel, fontsize='14')
    for i, v in enumerate(C_Type):
        plt.text(i*1.2 - 0.15, v+0.2, f'{v:.2f}',
                 color='blue', fontsize='11', fontweight='bold')
    plt.show()
    return Type_


# def plot_system_summary(system):
#     # Get data
#     df = report_cost(system.units)
#     Type = df['Type']
#     CI_unit = df['Capital Investment (10^6 USD)']
#     CU_unit = df['Utility Cost (10^6 USD/yr)']
#     CI_Type = []
#     CU_Type = []
#     t_old = Type[0]
#     r = 0
#     ci = 0
#     cu = 0
#     Type_ = []
#     for t in Type:
#         if t == t_old:
#             ci += CI_unit[r]
#             cu += CU_unit[r]
#         else:
#             t_old = t
#             CI_Type.append(ci)
#             CU_Type.append(cu)
#             Type_.append(t_old)
#             ci = CI_unit[r]
#             cu = CU_unit[r]
#         r += 1
    
#     # Make subplot
#     f, (ax1, ax2) = plt.subplots(2, 1, sharey=True)
#     pos = np.arange(len(Type_))
#     ax2.bar(pos, CI_Type, align='center', alpha=0.5)
#     plt.sca(ax2)
#     plt.xticks(pos, Type_)
#     plt.ylabel('Capital Investment (10^6 USD)')
#     plt.title('Cost Summary')
    
#     ax1.bar(pos, CU_Type, align='center', alpha=0.5)
#     plt.sca(ax1)
#     plt.ylabel('Capital Investment (10^6 USD)')
#     plt.title('Cost Summary')
#     plt.show()
    
#     return ax1, ax2