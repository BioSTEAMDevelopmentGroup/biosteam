# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst
try:
    from .profile import register
except:
    def register(*args, **kwargs):
        return lambda f: f

__all__ = (
    'create_system_butanol_purification',
)

@register(
    'butanol_purification', 'Butanol purification',
    2, [0.5, 1, 1.5, 2], 'BtOH\nsep.'
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
        x_bot=0.0001, y_top=0.2, k=1.2, Rmin=0.01,
        LHK=('Butanol', 'Water'),
    )
    
    # Decanter
    settler = bst.StageEquilibrium(
        'settler',
        ins=(distillate_0, distillate_1), 
        outs=('butanol_rich', water_rich),
        phases=('L', 'l'),
        top_chemical='Butanol',
        T=310,
    )
    
    # Butanol purification
    butanol_distiller = bst.BinaryDistillation(
        ins=settler-0,
        outs=(distillate_1, 'butanol'),
        x_bot=0.0001, y_top=0.6, k=1.2, Rmin=0.01,
        LHK=('Water', 'Butanol'),
    )
    
    sys = bst.System.from_units(units=[water_distiller, settler, butanol_distiller], algorithm=alg)
    return sys
