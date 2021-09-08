# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2021-2022, Yalin Li <zoe.yalin.li@gmail.com>, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

def test_qsdsan():
    from exposan import bwaise as bw
    # TODO: Make tests with assertions
    # Just make sure it can simulate for now
    for i in (bw.sysA, bw.sysB, bw.sysC): i.simulate() 

if __name__ == '__main__':
    test_qsdsan() 