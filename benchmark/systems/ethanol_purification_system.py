# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 13:38:11 2024

@author: cortespea
"""
import biosteam as bst

__all__ = (
    'create_system_ethanol_purification',
)

def create_system_ethanol_purification(alg):
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    distilled_beer = bst.Stream(
        'distilled_beer', 
        phase='g', T=386.16, P=212782, 
        Water=1403, Ethanol=545.7, units='kmol/hr'
    )
    ethanol = bst.Stream(
        'ethanol'
    )
    recycle = bst.Stream('recycle')
  
    # Mix ethanol Recycle (Set-up)
    
    # D303 = bst.BinaryDistillation(
    #     ins=(distilled_beer, recycle),
    #     x_bot=3.9106e-06, y_top=0.80805, k=1.2, Rmin=0.01,
    #     LHK=('Ethanol', 'Water'),
    #     P=212782,
    #     is_divided=True
    # )
    
    D303 = bst.MESHDistillation(
        ins=(distilled_beer, recycle),
        LHK=('Ethanol', 'Water'),
        P=212782,
        N_stages=20,
        feed_stages=(11,),
        reflux=2.42,
        boilup=0.4036,
        use_cache=True,
    )
    
    
    # Molecular sieve
    U301 = bst.Separator(
        ins=D303-0,
        outs=(recycle, ethanol),
        split=(2165.14/13356.04, 1280.06/1383.85),
        order=('Ethanol', 'Water'),
        T=115+273.15,
    )
    sys = bst.System.from_units(units=[D303, U301], algorithm=alg)
    return sys

def test_ethanol_purification_system():
    po = create_system_ethanol_purification('phenomena-oriented')
    po.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    sm = create_system_ethanol_purification('sequential modular')
    sm.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=20)
    po.simulate()
    sm.simulate()
    
if __name__ == "__main__":
    test_ethanol_purification_system()