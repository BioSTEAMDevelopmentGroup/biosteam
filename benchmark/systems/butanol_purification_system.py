# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst

__all__ = (
    'create_system_butanol_purification',
)

def create_system_butanol_purification(alg):
    bst.settings.set_thermo(['Water', 'Butanol'], cache=True)
    feed = bst.Stream(
        'feed', 
        phase='l', 
        Water=900,
        Butanol=100,
        units='kmol/hr',
    )
    water_rich = bst.Stream('water_rich')
    distillate_0 = bst.Stream('distillate_0')
    distillate_1 = bst.Stream('distillate_1')
    
    # Dewater
    water_distiller = bst.BinaryDistillation(
        ins=(feed, water_rich), outs=(distillate_0, 'water'),
        # Hr=0.66, Lr=0.999,
        x_bot=0.0001, y_top=0.2, 
        k=1.2, Rmin=0.01,
        LHK=('Butanol', 'Water'),
        partial_condenser=False,
    )
    # water_distiller.decoupled = True
    
    # Decanter
    mixer = bst.Mixer('mixer', ins=(distillate_0, distillate_1))
    settler = bst.StageEquilibrium(
        'settler',
        ins=mixer-0, 
        outs=('butanol_rich', water_rich),
        phases=('L', 'l'),
        top_chemical='Butanol',
        T=310,
    )
    
    # Butanol purification
    butanol_distiller = bst.BinaryDistillation(
        ins=settler-0,
        outs=(distillate_1, 'butanol'),
        # Hr=0.61, Lr=0.999,
        x_bot=0.0001, y_top=0.6, 
        k=1.2, Rmin=0.01,
        LHK=('Water', 'Butanol'),
        partial_condenser=False,
    )
    # butanol_distiller.decoupled = True
    
    sys = bst.System.from_units(units=[water_distiller, mixer, settler, butanol_distiller], algorithm=alg)
    # breakpoint()
    return sys

def test_butanol_purification_system():
    po = create_system_butanol_purification('phenomena-oriented')
    po.flatten()
    po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_system_butanol_purification('sequential modular')
    sm.flatten()
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    po.simulate()