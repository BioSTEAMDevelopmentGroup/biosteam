# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2021-2022, Yalin Li <zoe.yalin.li@gmail.com>,
#                          Yoel Cortes-Pena <yoelcortes@gmail.com>,
#                          Joy Zhang <joycheung1994@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

def test_dyn_sys():
    from qsdsan import processes as pc, sanunits as su, set_thermo, System
    import numpy as np
    from numpy.testing import assert_allclose
    
    cmps = pc.load_asm1_cmps()
    set_thermo(cmps)
    DI = su.DynamicInfluent('Dyn_Inf')
    S1 = su.Splitter('Split', ins=DI-0, split=0.3, init_with='WasteStream')
    M1 = su.Mixer('Mix', ins=(S1-0, S1-1), outs=('Dyn_Eff'))
    sys = System('test_sys', path=(DI, S1, M1))
    sys.set_dynamic_tracker(DI.outs[0], M1.outs[0])

    t = 1
    t_step = 0.05
    sys.simulate(t_span=(0,t),
                 t_eval=np.arange(0, t+t_step, t_step))
    dinf = sys.units[0].outs[0]
    deff = sys.units[-1].outs[0]
    assert_allclose(deff.scope.record, dinf.scope.record, rtol=1e-12)


def test_qsdsan():
    from exposan import bwaise as bw
    for i in (bw.sysA, bw.sysB, bw.sysC): i.simulate() 
    test_dyn_sys()

if __name__ == '__main__':
    test_qsdsan() 