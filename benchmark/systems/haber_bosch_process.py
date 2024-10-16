# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import thermosteam as tmo
import biosteam as bst
import numpy as np
from thermosteam.constants import R
from math import exp

__all__ = (
    'create_system_haber_bosch_process',
)

def create_system_haber_bosch_process(alg='sequential modular'):
    bst.settings.set_thermo(['N2', 'H2', 'NH3'], cache=True)
    bst.settings.mixture.include_excess_energies = True
    
    with bst.System(algorithm=alg) as sys:
        feed = bst.Stream(
            'feed',
            H2=3,
            N2=1,
            P=300e5,
            phase='g'
        )
        preheater = bst.HXutility(ins=feed, T=450 + 273.15)
        reactor = bst.ReactivePhaseStage(
            ins=preheater-0, T=450 + 273.15, P=300e5,
            reaction=bst.Reaction('N2 + H2 -> NH4', reactant='N2', X=0.15, correct_atomic_balance=True),
        )
        
        makeup_butanol = bst.Stream('makeup_butanol', Butanol=1)
        makeup_butanol.T = makeup_butanol.bubble_point_at_P().T
        recycle_butanol = bst.Stream('recycle_butanol')
        esterification_reflux = bst.Stream('esterification_reflux')
        esterification = bst.MESHDistillation(
            'esterification',
            ins=(feed, makeup_butanol, recycle_butanol, esterification_reflux), 
            outs=('empty', 'bottoms', 'esterification_distillate'),
            N_stages=6,
            feed_stages=(3, 3, 3, 0),
            stage_specifications={
                5: ('Boilup', 1),
                0: ('Boilup', 0),
            },
            liquid_side_draws={
                0: 1,
            },
            stage_reactions={
                i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
                for i in range(1, 5)
            },
            maxiter=10,
            # LHK=('LacticAcid', 'Butanol'),
            LHK=('Butanol', 'ButylLactate'),
            P=0.3 * 101325,
        )
        @esterification.add_specification(run=True)
        def adjust_flow():
            target = 5.85
            makeup_butanol.imol['Butanol'] = max(target - recycle_butanol.imol['Butanol'], 0)
        
        esterification.simulate()
        for i in esterification.stages: print(i.Hnet) 
        esterification.show()
        breakpoint()
        # esterification.simulate()
        # for i in esterification.stages: print(i.Hnet) 
        # esterification.stage_reactions={
        #         i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
        #         for i in range(1, 17)
        #     }
        # esterification.LHK=('Butanol', 'ButylLactate')
        # breakpoint()
        # esterification.simulate()
        esterification_settler = bst.StageEquilibrium(
            'esterification_settler',
            ins=(esterification-2), 
            outs=(esterification_reflux, 'water_rich'),
            phases=('L', 'l'),
            top_chemical='Butanol',
        )
        water_distiller = bst.BinaryDistillation(
            ins=esterification_settler-1, outs=('water_rich_azeotrope', 'water'),
            x_bot=0.0001, y_top=0.2, k=1.2, Rmin=0.01,
            LHK=('Butanol', 'Water'),
        )
        splitter = bst.Splitter(ins=water_distiller-1, split=0.5) # TODO: optimize split
        hydrolysis_reflux = bst.Stream('hydrolysis_reflux')
        hydrolysis = bst.MESHDistillation(
            'hydrolysis',
            ins=(esterification-1, splitter-0, hydrolysis_reflux),
            outs=('empty', 'lactic_acid', 'hydrolysis_distillate'),
            N_stages=53,
            feed_stages=(27, 50, 0),
            stage_specifications={
                0: ('Boilup', 0),
                52: ('Boilup', 1),
            },
            liquid_side_draws={
                0: 1.0,
            },
            stage_reactions={
                i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
                for i in range(1, 52) # It will run in reverse
            },
            P=101325,
            LHK=('Butanol', 'LacticAcid'),
        )
        
        # @esterification.add_specification(run=True)
        # def adjust_flow():
        #     target = 5.85
        #     makeup_butanol.imol['Butanol'] = max(target - recycle_butanol.imol['Butanol'], 0)
        
        # Decanter
        butanol_rich_azeotrope = bst.Stream('butanol_rich_azeotrope')
        hydrolysis_settler = bst.StageEquilibrium(
            'settler',
            ins=(hydrolysis-2, water_distiller-0, butanol_rich_azeotrope), 
            outs=('butanol_rich_extract', hydrolysis_reflux),
            phases=('L', 'l'),
            top_chemical='Butanol',
            T=310,
        )
        
        # Butanol purification
        butanol_distiller = bst.BinaryDistillation(
            ins=(hydrolysis_settler-0),
            outs=(butanol_rich_azeotrope, recycle_butanol),
            x_bot=0.0001, y_top=0.6, k=1.2, Rmin=0.01,
            LHK=('Water', 'Butanol'),
        )
        
    return sys

if __name__ == '__main__':
    sys = create_system_acetic_acid_reactive_purification()
    sys.flatten()
    sys.diagram()
    sys.simulate()