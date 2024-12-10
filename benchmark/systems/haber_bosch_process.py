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
    chemicals = bst.Chemicals(['N2', 'H2', 'NH3'], cache=True)
    bst.settings.set_thermo(
        chemicals, 
        # mixture=bst.PRMixture.from_chemicals(chemicals),
        # Phi=bst.PRFugacityCoefficients,
        # Gamma=bst.PRActivityCoefficients,
    )
    bst.settings.mixture.include_excess_energies = True
    
    with bst.System(algorithm=alg) as sys:
        feed = bst.Stream(
            'feed',
            H2=3,
            N2=1,
            P=300e5,
            phase='g'
        )
        recycle = bst.Stream()
        P_rxn = 300e5
        P_condenser_drop = 2e5
        P_preheater_drop = 2e5
        preheater = bst.SinglePhaseStage(ins=[feed, recycle], T=450 + 273.15, P=P_rxn)
        @preheater.add_specification(run=True)
        def adjust_fresh_feed():
            feed.imol['N2'] = 1 - recycle.imol['N2']
            feed.imol['H2'] = 3 - recycle.imol['H2']
            # if sys.algorithm[0].lower() == 'p': breakpoint()
            
        @feed.material_balance
        def fresh_feed_flow_rate():
            f = np.ones(chemicals.size)
            r = np.zeros(chemicals.size)
            v = np.zeros(chemicals.size)
            index = chemicals.indices(['N2', 'H2'])
            v[index] = [1, 3]
            r[index] = [1, 1]
            return (
                {feed: f,
                 recycle: r},
                 v
            )
            
        reactor = bst.ReactivePhaseStage(
            ins=preheater-0, T=450 + 273.15, P=P_rxn,
            reaction=bst.Reaction('N2 + H2 -> NH3', reactant='N2', X=0.15, correct_atomic_balance=True),
        )
        flash = bst.StageEquilibrium(ins=reactor-0, T=273.15 - 20, P=P_rxn - P_condenser_drop, phases=('g', 'l'))
        compressor = bst.IsentropicCompressor(ins=flash-0, outs=recycle, P=P_rxn + P_preheater_drop)
        
    return sys

if __name__ == '__main__':
    sm = create_system_haber_bosch_process()
    sm.diagram()
    sm.set_tolerance(mol=1e-9, rmol=1e-9, method='fixed-point')
    sm.simulate()
    
    print('-----------')
    po = create_system_haber_bosch_process('Phenomena-oriented')
    po.diagram()
    po.set_tolerance(mol=1e-9, rmol=1e-9, method='fixed-point')
    po.simulate()
    