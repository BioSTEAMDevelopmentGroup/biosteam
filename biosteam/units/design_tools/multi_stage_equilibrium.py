# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
General functional algorithms based on MESH equations to solve multi-stage
equilibrium.

"""
import numpy as np

def single_component_flow_rates_for_multi_stage_extration_without_side_draws(
        N_stages,
        phase_ratios,
        partition_coefficients, 
        feed, 
        solvent
    ):
    """
    Solve flow rates for a single component across a multi stage extraction
    operation without side draws. 

    Parameters
    ----------
    N_stages : int
        Number of stages.
    phase_ratios : 1d array
        Phase ratios by stage. The phase ratio for a given stage is 
        defined as F_a / F_b; where F_a and F_b are the flow rates 
        of phase a and b leaving the stage respectively.
    partition_coefficients : 1d array
        Partition coefficients by stage. The partition coefficient for a given
        stage is defined as x_a / x_b; where x_a and x_b are the fraction of
        the component in phase a and b leaving the stage.
    feed : float
        Component flow rate in feed entering stage 1.
    solvent : float
        Component flow rate in solvent entering stage N.

    Returns
    -------
    phase_a_flow_rates : 1d array
        Component flow rates in phase a by stage.
    phase_b_flow_rates : 1d array
        Component flow rates in phase b by stage.

    """
    component_ratios = phase_ratios * partition_coefficients
    A = np.eye(N_stages) * (1 + component_ratios) 
    stages_index = np.arange(N_stages)
    A[stages_index[:-1], stages_index[1:]] = -1
    A[stages_index[1:], stages_index[:-1]] = -component_ratios[:-1]
    b = np.zeros(N_stages)
    b[0] = feed
    b[N_stages] = solvent
    phase_b_flow_rates = np.solve(A, b)
    phase_a_flow_rates = phase_ratios * phase_b_flow_rates
    return phase_a_flow_rates, phase_b_flow_rates

def flow_rates_for_multi_stage_extration_without_side_draws(
        N_stages,
        phase_fractions,
        partition_coefficients, 
        feed, 
        solvent
    ):
    """
    Solve flow rates for a single component across a multi stage extraction
    operation without side draws. 

    Parameters
    ----------
    N_stages : int
        Number of stages.
    phase_fractions : 1d array
        Phase fractions by stage. The phase fraction for a given stage is 
        defined as F_a / (F_a + F_b); where F_a and F_b are the flow rates 
        of phase a and b leaving the stage respectively.
    partition_coefficients : Iterable[1d array]
        Partition coefficients with stages by row and components by column.
        The partition coefficient for a component in a given stage is defined 
        as x_a / x_b; where x_a and x_b are the fraction of the component in 
        phase a and b leaving the stage.
    feed : Iterable[float]
        Flow rates of all components in feed entering stage 1.
    solvent : Iterable[float]
        Flow rates of all components in solvent entering stage N.

    Returns
    -------
    phase_a_flow_rates : tuple[1d array]
        Flow rates in phase a with stages by row and components by column.
    phase_b_flow_rates : tuple[1d array]
        Flow rates in phase b with stages by row and components by column.

    """
    phase_ratios = phase_fractions / (1. - phase_fractions)
    phase_a_flow_rates, phase_b_flow_rates = zip(*[
        single_component_flow_rates_for_multi_stage_extration_without_side_draws(
            N_stages, phase_ratios, K, f, s) 
        for K, f, s in zip(partition_coefficients, feed, solvent)
    ])
    return phase_a_flow_rates, phase_b_flow_rates
