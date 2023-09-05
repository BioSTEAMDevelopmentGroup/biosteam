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
    feed = bst.Stream('feed', Water=100, Ethanol=100)
    recycle = bst.Stream('recycle', Water=100, Ethanol=100)
    M1 = bst.Mixer(ins=[feed, recycle])
    F1 = bst.Flash(ins=M1.outlet, outs=['vapor_product', 'liquid'], T=365, P=101325)
    S1 = bst.Splitter(ins=F1.liquid, outs=[recycle, 'liquid_product'], split=0.4)
    
    sys = bst.System.from_units('sys', [M1, F1, S1])
    sys.set_tolerance(mol=1e-6, rmol=1e-6, rT=1e-6, T=1e-6)
    model = bst.Model(sys)
    @model.parameter(distribution=shape.Uniform(0.1, 0.9), kind='coupled')
    def set_ethanol_fraction(x_ethanol):
        total_flow = feed.F_mol
        feed.imol['Ethanol'] = total_flow * x_ethanol
        feed.F_mol = total_flow
        
    convergence_model = bst.NullConvergenceModel(
        predictors=[set_ethanol_fraction], 
    )
    model.load_samples(model.sample(100, rule='L', seed=1), optimize=True)
    model.evaluate(design_and_cost=False, convergence_model=convergence_model)
    R2_null, _ = convergence_model.R2()
        
    convergence_model = bst.ConvergenceModel(
        predictors=[set_ethanol_fraction], local_weighted=False,
        model_type=bst.InterceptLinearRegressor, save_prediction=True,
    )
    model.evaluate(design_and_cost=False, convergence_model=convergence_model)
    summary, _ = convergence_model.R2()
    R2p = summary['predicted']
    R2f = summary['fitted']
    assert R2f['min'] > R2p['min'] > R2_null['min']
    assert R2f['max'] > R2p['max'] > R2_null['max']
    
if __name__ == '__main__':
    test_convergence_model()
    