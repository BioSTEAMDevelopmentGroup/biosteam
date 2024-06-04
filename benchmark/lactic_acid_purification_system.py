# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst
from thermosteam.constants import R
from math import exp
try:
    from .profile import register
except:
    def register(*args, **kwargs):
        return lambda f: f

__all__ = (
    'create_system_lactic_acid_purification',
)

@register(
    'lactic_acid_purification', 'Lactic acid purification',
    10, [2, 4, 6, 8, 10], 'LA\nsep.'
)
def create_system_lactic_acid_purification(alg='sequential modular'):
    bst.settings.set_thermo(['Water', 'LacticAcid', 'ButylLactate', 'Butanol'], cache=True)
    
    class Esterification(bst.KineticReaction):
        
        def volume(self, stream): # kg of catalyst
            rho_cat = 770 # kg / m3
            liquid_volume = self.liquid_volume
            catalyst_volume = 0.5 * liquid_volume
            catalyst_mass = catalyst_volume * rho_cat
            return catalyst_mass
        
        def rate(self, stream): # kmol/kg-catalyst/hr
            T = stream.T
            # if T > 370: return 0 # Prevents multiple steady states
            kf = 2.59e4 * exp(-5.340e4 / (R * T))
            kr = 3.80e3 * exp(-5.224e4 / (R * T))
            H2O, LA, BuLA, BuOH = stream.mol / stream.F_mol
            return self.stoichiometry * 3600 * (kf * LA * BuOH - kr * BuLA * H2O) # kmol / kg-catalyst / hr
    
    with bst.System(algorithm=alg) as sys:
        feed = bst.Stream(
            'feed',
            LacticAcid=4.174,
            Water=5.470,
        )
        makeup_butanol = bst.Stream('makeup_butanol')
        recycle_butanol = bst.Stream('recycle_butanol')
        esterification_reflux = bst.Stream('esterification_reflux')
        esterification = bst.MESHDistillation(
            'esterification',
            ins=(feed, makeup_butanol, recycle_butanol, esterification_reflux), 
            outs=('empty', 'bottoms', 'esterification_distillate'),
            N_stages=17,
            feed_stages=(1, 16, 16, 0),
            stage_specifications={
                16: ('Boilup', 1),
                0: ('Boilup', 0),
            },
            liquid_side_draws={
                0: 1.0,
            },
            stage_reactions={
                i: Esterification('LacticAcid + Butanol -> Water + ButylLactate', reactant='LacticAcid')
                for i in range(1, 17)
            },
            maxiter=200,
            LHK=('Butanol', 'ButylLactate'),
            P=0.3 * 101325,
        )
        @esterification.add_specification(run=True)
        def adjust_flow():
            target = 5.85
            makeup_butanol.imol['Butanol'] = max(target - recycle_butanol.imol['Butanol'], 0)
            
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