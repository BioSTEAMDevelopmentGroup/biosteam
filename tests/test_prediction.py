# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Yalin Li <zoe.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

def test_convergence_model():
    import biosteam as bst
    from chaospy import distributions as shape
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream(Water=100, Ethanol=100)
    recycle = bst.Stream(Water=100, Ethanol=100)
    
    with bst.System() as sys:
        M1 = bst.Mixer(ins=[feed, recycle])
        F1 = bst.Flash(ins=M1.outlet, outs=['vapor', 'liquid_product'], V=0.5, P=101325)
        S1 = bst.Splitter(ins=F1.vapor, outs=['vapor_product', recycle], split=0.5)
    
    sys.set_tolerance(mol=1e-6, rmol=1e-6, rT=1e-6, T=1e-6, method='aitken', maxiter=100)
    model = bst.Model(sys)
    total_flow = feed.F_mol
    @model.parameter(distribution=shape.Uniform(0.1, 0.4), kind='coupled')
    def set_ethanol_fraction(x_ethanol):
        feed.imol['Ethanol'] = total_flow * x_ethanol
        feed.imol['Water'] = total_flow * (1 - x_ethanol)
        
    # Compare linear model against unsorted
    convergence_model = bst.NullConvergenceModel(
        parameters=[set_ethanol_fraction], 
        save_prediction=True,
        system=sys,
    )
    model.load_samples(model.sample(100, rule='L', seed=1), sort=False)
    model.evaluate(design_and_cost=False, convergence_model=convergence_model)
    R2_null, _ = convergence_model.R2()
        
    convergence_model = bst.ConvergenceModel(
        parameters=[set_ethanol_fraction], local_weighted=False,
        model_type=bst.InterceptLinearRegressor, save_prediction=True,
        system=sys, 
    )
    model.evaluate(design_and_cost=False, convergence_model=convergence_model)
    summary, _ = convergence_model.R2()
    R2p = summary['predicted']
    R2f = summary['fitted']
    for key in ('min', 'mean', 'max'):
        assert R2p[key] > R2_null[key]
        assert R2f[key] > R2_null[key]
    

if __name__ == '__main__':
    test_convergence_model()
    