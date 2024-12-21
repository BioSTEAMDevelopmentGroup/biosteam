# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
from numpy.testing import assert_allclose


def test_fit():
    import biosteam as bst
    import numpy as np
    Ce = np.array([
        0.06013, 0.07133, 0.07133, 0.31769, 0.31769, 0.31769,
        0.63124, 0.64244, 0.63124, 2.03102, 2.03102, 2.03102,
        4.16988, 4.16988, 4.16988, 5.00974, 4.99854, 4.99854
    ])
    qe = np.array([
        1.24, 1.24, 1.24, 1.87, 1.87, 1.87, 2.85, 2.84, 2.85,
        4.48, 4.48, 4.48, 5.5, 5.5, 5.5, 6.53, 6.56, 6.56
    ])
    t = np.array([0.0, 31.0, 60.0, 90.0, 120.0, 150.0]) / 60
    C = np.array([7.0, 2.45, 1.03, 0.43, 0.26, 0.13])
    volume = 0.075
    adsorbent = 1
    bst.AdsorptionColumn.plot_isotherm_and_mass_transfer_coefficient_fit(
        Ce, qe, t, C, volume, adsorbent
    )
    
def test_adsorption():
    import biosteam as bst
    chemicals = bst.Chemicals([
        'Xylene', 
        bst.Chemical('Ink', search_db=False, default=True, phase='l'),
        bst.Chemical('ActivatedCarbon', search_db=False, default=True, phase='l')
    ])
    bst.settings.set_thermo(chemicals)
    chemicals.compile()
    chemicals.define_group('Solvent', ['Xylene'])
    feed = bst.Stream(ID='feed', phase='l', T=298, P=1.01e+06,
                      Xylene=1000, Ink=0.001, units='kg/hr')
    # feed.T = 450
    # Adsorption equilibrium and rate data from Dr. Tianwei Yan from Prof. George Huber's group.
    A1 = bst.AdsorptionColumn(
        'A1', ins=feed, outs='effluent',
        cycle_time=1000,
        superficial_velocity=4,
        # isotherm_model='Freundlich',
        # isotherm_args = (3.31e3, 2.58), 
        # k = 0.156, # [1 / hr]
        # Equilibrium constant [L/g] for Langmuir isotherm (multiply by 1000 from L/mg)
        # Maximum equilibrium loading [mg/g] for Langmuir isotherm
        isotherm_model='Langmuir',
        isotherm_args = (1.27e3, 6.99), 
        k = 0.3, #[1 / hr]
        void_fraction=1 - 0.38 / 0.8,  # Solid by wt [%]
        C_final_scaled=0.05,
        adsorbate='Ink',
    )
    A1.simulate()
    assert_allclose(A1.outs[0].mol, [9.419, 0., 0.], rtol=1e-3)
    assert 1e-3 < A1.CL_scaled[-1, -1] < 2 * A1.C_final_scaled # Shape of MTZ is not exact
    assert_allclose(A1.CL_scaled[-1, -1], 0.0842941306073549, rtol=1e-3)
    assert_allclose(A1.LES, 2.504118270726701)
    assert_allclose(A1.q0, 3.6870766727528457)
    assert_allclose(A1.design_results['Diameter'], 1.9764225277140537)
    assert_allclose(A1.get_design_result('Length', 'm'), A1.LES + A1.LUB)

def test_adsorption_bed_pressure_drop():
    f = lambda L, u: bst.unuts.adsorption.adsorption_bed_pressure_drop(u=u/3600, L=L)
    X, Y, Z = bst.plots.generate_contour_data(f, [1, 10], [4, 14], n=10)
    


if __name__ == '__main__':
    test_fit()
    test_adsorption_bed_pressure_drop()
