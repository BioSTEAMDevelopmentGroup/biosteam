# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 08:32:10 2025

@author: yoelr
"""
import biosteam as bst

__all__ = (
    'create_stripper_system', 
    'create_flash_system',
)

def create_stripper_system(alg):
    bst.settings.set_thermo(['AceticAcid', 'Water', 'MTBE'], cache=True)
    feed = bst.Stream('feed', Water=75, AceticAcid=5, MTBE=20, T=320)
    steam = bst.Stream('steam', Water=100, phase='g', T=390)
    stripper = bst.Stripper('D1',
        N_stages=5, ins=[feed, steam], 
        outs=['vapor', 'liquid'],
        solute="AceticAcid", 
    )
    return bst.System.from_units(units=[stripper], algorithm=alg)
    
def create_flash_system(alg):
    bst.settings.set_thermo(
        ['Water', 'Ethanol'], cache=True
    )
    with bst.System(algorithm=alg) as system:
        bst.StageEquilibrium(
            ins=bst.Stream('feed', Water=1, Ethanol=1), 
            outs=['vapor', 'liquid'],
            phases=('g', 'l')
        )
    return system