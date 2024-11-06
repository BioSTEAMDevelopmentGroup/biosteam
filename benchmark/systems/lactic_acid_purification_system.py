# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst
from thermosteam.constants import R
from math import exp

__all__ = (
    'create_system_lactic_acid_purification',
)

def create_system_lactic_acid_purification(alg='sequential modular'):
    bst.settings.set_thermo(['Water', 'LacticAcid', 'EthylLactate', 'Ethanol', 'SuccinicAcid'], cache=True)
    
    class Esterification(bst.KineticReaction):
        catalyst_fraction = 0.5 
        
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = self.catalyst_fraction * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream):
            T = stream.T
            if T > 365: return 0 # Prevents multiple steady states.
            kf = 6.52e3 * exp(-4.8e4 / (R * T))
            kr = 2.72e3 * exp(-4.8e4 / (R * T))
            LaEt, La, H2O, EtOH, _ = stream.mol / stream.F_mol
            return 3600 * (kf * La * EtOH - kr * LaEt * H2O) # kmol / kg-catalyst / hr
    
    
    with bst.System(algorithm=alg) as sys:
        feed = bst.Stream(
            'feed',
            LacticAcid=4.174,
            Water=5.460,
            SuccinicAcid=0.531,
            EthylLactate=1e-6,
            total_flow=10.165,
            P=0.106 * 101325,
            T=72.5 + 273.15
        )
        feed.T = feed.bubble_point_at_P(P=0.2 * 101325).T
        makeup_ethanol = bst.Stream('makeup_ethanol', Ethanol=0.035, P=2 * 101325, phase='g')
        makeup_ethanol.T = makeup_ethanol.dew_point_at_P(P=0.2 * 101325).T
        recycle_ethanol = bst.Stream('recycle_ethanol')
        esterification = bst.MESHDistillation(
            'esterification',
            ins=(feed, makeup_ethanol, recycle_ethanol), 
            outs=('empty', 'bottoms', 'esterification_disti2llate'),
            N_stages=29,
            feed_stages=(1, 27, 27),
            stage_specifications={
                0: ('Reflux', 1),
                28: ('Boilup', 1),
                # 0: ('Temperature', 46.8 + 273.15),
                # 28: ('Temperature', 229.6 + 273.15),
            },
            full_condenser=True,
            stage_reactions={
                i: Esterification('LacticAcid + Ethanol -> Water + EthylLactate', reactant='LacticAcid')
                for i in range(1, 28)
            },
            maxiter=50,
            LHK=('Ethanol', 'EthylLactate'),
            P=0.2 * 101325,
            use_cache=True
        )
        @esterification.add_specification(run=True)
        def adjust_flow():
            target = 0.035 + 16.653
            makeup_ethanol.imol['Ethanol'] = max(target - recycle_ethanol.imol['Ethanol'], 0)
        bst.PhasePartition.B_relaxation_factor = 0.9
        bst.PhasePartition.K_relaxation_factor = 0
        catalyst_fraction = 0
        dc = 1e-6
        while catalyst_fraction < 0.5:
            dc *= 1.5
            if dc > 5e-3: dc = 1e-3
            if dc < 1e-4: dc = 1e-4
            catalyst_fraction += dc
            if catalyst_fraction > 0.5:
                catalyst_fraction = 0.5
            print('----')
            print(catalyst_fraction)
            print('----')
            Esterification.catalyst_fraction = catalyst_fraction
            esterification.simulate()
            for i in esterification.stages: print(i.Hnet) 
            esterification.show()
        
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
    sys = create_system_lactic_acid_purification()
    sys.flatten()
    sys.diagram()
    sys.simulate()