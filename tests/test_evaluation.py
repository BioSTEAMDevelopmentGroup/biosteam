# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>, Yalin Li <zoe.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
from numpy.testing import assert_allclose
import numpy as np

def test_parameter_hook():
    import biosteam as bst
    from chaospy.distributions import Uniform
    sys = bst.System(None, ())
    model = bst.Model(sys)
    
    @model.parameter(distribution=Uniform(0., 1.), hook=round)
    def set_param(x):
        pass
    
    np.random.seed(0)
    samples = model.sample(10, 'L')
    model.load_samples(samples)
    actual_values = [[1.],
                     [0.],
                     [0.],
                     [0.],
                     [0.],
                     [0.],
                     [1.],
                     [1.],
                     [1.],
                     [1.]]
    assert_allclose(model.table.values, actual_values) 

def create_evaluation_model(cache=[None]):
    if cache[0] is not None: return cache[0]
    import biosteam as bst
    from chaospy.distributions import Uniform
    sys = bst.System(None, ())
    model = bst.Model(sys)
    
    parameter_box = [100.]
    @model.parameter(bounds=[90, 110],
                     distribution=Uniform(90, 110),
                     units='kg/hr')
    def set_parameter(parameter_1):
        parameter_box[0] = parameter_1
    
    @model.parameter(bounds=[0., 1.],
                     distribution=Uniform(0., 1.))
    def set_parameter(parameter_2): ...
    
    @model.metric
    def good_metric():
        p = parameter_box[0]
        return 2. * p * p / (p + 100.) 
    
    switch_box = [False]
    @model.metric
    def bad_metric():
        give_nan = switch_box[0]
        switch_box[0] = not give_nan
        p = parameter_box[0]
        return float('Nan') if give_nan else 2. * p * p / (p + 100.) 
    
    @model.metric
    def non_correlated_metric():
        return float(switch_box[0])
    
    sample = model.sample(100, 'L')
    model.load_samples(sample)
    model.evaluate()
    cache[0] = model
    return model
    
def test_pearson_r():
    model = create_evaluation_model()
    with pytest.raises(ValueError):
        rho, p = model.pearson_r(filter='none')
    
    rho, p = model.pearson_r(filter='propagate nan')
    NaN = float('NaN')
    expected = np.array([[1., NaN, 0.],
                         [0,  NaN, 0.]])
    index = ~np.isnan(expected)
    assert_allclose(np.round(rho.values[index]), expected[index], atol=0.15)
    
    with pytest.raises(ValueError):
        rho, p = model.pearson_r(filter='raise nan error')
    
    rho, p = model.pearson_r(filter='omit nan')
    expected = np.array([[1., 1., 0.],
                         [0,  0., 0.]])
    assert_allclose(np.round(rho), expected, atol=0.15)

def test_spearman_r():
    model = create_evaluation_model()
    rho, p = model.spearman_r()
    NaN = float('NaN')
    expected = np.array([[1., NaN, 0.],
                         [0,  NaN, 0.]])
    index = ~np.isnan(expected)
    assert_allclose(rho.values[index], expected[index], atol=0.2)
    
    rho, p = model.spearman_r(filter='omit nan')
    expected = np.array([[1., 1., 0.],
                         [0,  0., 0.]])
    assert_allclose(np.round(rho), expected, atol=0.15)
    
def test_kendall_tau():
    model = create_evaluation_model()
    rho, p = model.kendall_tau()
    NaN = float('NaN')
    expected = np.array([[1., NaN, 0.],
                         [0,  NaN, 0.]])
    index = ~np.isnan(expected)
    assert_allclose(rho.values[index], expected[index], atol=0.15)
    
    tau, p = model.kendall_tau(filter='omit nan')
    expected = np.array([[1., 1., 0.],
                         [0,  0., 0.]])
    assert_allclose(np.round(tau), expected, atol=0.15)
    
def test_model_index():
    import biosteam as bst
    bst.settings.set_thermo(['Water'])
    with bst.System() as sys:
        H1 = bst.HXutility('H1', ins=bst.Stream('feed', Water=1000), T=310)
    
    # Make sure repeated metrics raise an error
    biorefinery = bst.process_tools.UnitGroup('Biorefinery', sys.units)
    heating_duty = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    heating_duty_repeated = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    with pytest.raises(ValueError):
        model = bst.Model(sys, [heating_duty, heating_duty_repeated])
    
    model = bst.Model(sys, [heating_duty])
    
    # Make sure repeated parameters raise an error
    @model.parameter(element=H1)
    def set_temperature(temperature):
        H1.T = temperature
    
    with pytest.raises(ValueError):
        @model.parameter(element=H1)
        def set_efficiency(temperature):
            H1.T = temperature
    
    bst.default()
    
def test_model_sample():
    model = create_evaluation_model()
    # Just make sure they run
    for rule in ('MORRIS', 'FAST', 'RBD', 'SOBOL'): model.sample(100, rule)
    problem = model.problem()
    for rule in ('MORRIS', 'FAST', 'RBD', 'SOBOL'): model.sample(100, rule, problem=problem)

def test_copy():
    model = create_evaluation_model()
    copy = model.copy()
    assert copy.table is not None
    assert copy.table is not model.table
    assert copy.table.shape == model.table.shape
    copy.evaluate() # Make sure it runs
    
def test_model_exception_hook():
    import biosteam as bst
    import pytest
    from biorefineries import sugarcane as sc
    from chaospy import distributions as shape
    from warnings import simplefilter
    import numpy as np
    bst.settings.set_thermo(sc.chemicals)
    simplefilter("ignore")
    IRR_metric = bst.Metric('Internal rate of return', sc.sugarcane_tea.solve_IRR)
    metrics = [IRR_metric]
    sugarcane_model = bst.Model(sc.sugarcane_sys, metrics)
    baseline = sc.sugarcane.F_mass
    distribution = shape.Triangle(-baseline , baseline , 2*baseline) # Negative value should fail
    
    @sugarcane_model.parameter(element=sc.sugarcane, distribution=distribution, units='kg/hr')
    def set_sugarcane_flow_rate(flow_rate):
        sc.sugarcane.F_mass = flow_rate
    
    np.random.seed(0)
    samples = sugarcane_model.sample(5, 'L')
    sugarcane_model.load_samples(samples)
    
    # Without an exception hook, the same behavior will result (NaN values for failed evaluations)
    sugarcane_model.evaluate(convergence_model=False)
    assert np.isnan(sugarcane_model.table.values).any()
    
    InfeasibleRegion = bst.exceptions.InfeasibleRegion
    
    # This will provide a more understandable IRR result for infeasible regions
    def exception_hook(exception, sample): 
        if isinstance(exception, (InfeasibleRegion, ValueError, RuntimeError)):
            return [0]
        else: raise exception
    sugarcane_model.exception_hook = exception_hook
    sugarcane_model.evaluate(convergence_model=False)
    assert not np.isnan(sugarcane_model.table.values).any()
    
    # This will raise an exception due to negative flow rates
    def exception_hook(exception, sample): 
        if isinstance(exception, InfeasibleRegion):
            raise exception
    sugarcane_model.exception_hook = exception_hook
    with pytest.raises(InfeasibleRegion): sugarcane_model.evaluate(convergence_model=False)
    
    # This will raise an exception regardless
    def exception_hook(exception, sample): 
        raise exception
    sugarcane_model.exception_hook = exception_hook
    with pytest.raises(InfeasibleRegion): sugarcane_model.evaluate(convergence_model=False)
    
    # Here is another cool thing we could do in the case where 
    # some metrics are expected to fail
    bad_metric = bst.Metric('bad metric', lambda: 1/0)
    sugarcane_model.metrics = (IRR_metric, bad_metric)
    sugarcane_model.load_samples(samples) # Metrics changed, so need to reload sample
    def exception_hook(exception, sample): 
        if not isinstance(exception, ZeroDivisionError): return
        sc.sugarcane_sys.simulate()
        values = []
        for i in sugarcane_model.metrics:
            try: x = i()
            except: x = None
            values.append(x)
        return values
    sugarcane_model.exception_hook = exception_hook
    sugarcane_model.evaluate(convergence_model=False)
    bad_metric_results = sugarcane_model.table[bad_metric.index]
    IRR_metric_results = sugarcane_model.table[IRR_metric.index]
    assert np.isnan(bad_metric_results).all()
    assert not np.isnan(IRR_metric_results).all()
    
    # Note that the possibilities are infinite here...
    
def test_kolmogorov_smirnov_d():
    # Test made by Yalin.
    import numpy as np
    import biosteam as bst
    from chaospy import distributions as shape
    
    bst.settings.set_thermo(['H2O', 'Ethanol'], cache=True)
    
    s1 = bst.Stream('s1', H2O=100)
    M1 = bst.MixTank(ins=s1)
    M2 = bst.MixTank(ins=M1-0)
    sys = bst.System('sys', path=(M1, M2))
    
    model = bst.Model(sys)
    
    baseline = 1
    distribution = shape.Uniform(lower=0.5, upper=1.5)
    @model.parameter(name='M1 tau', element=M1, kind='coupled', distribution=distribution,
                     units='hr', baseline=baseline)
    def set_M1_tau(i):
        M1.tau = i
        
    baseline = 1.75
    distribution = shape.Uniform(lower=1, upper=2)
    @model.parameter(name='M2 tau', element=M2, kind='coupled', distribution=distribution,
                     units='hr', baseline=baseline)
    def set_M2_tau(i):
        M2.tau = i
    
    model.metrics = [
        bst.Metric(name='tau1', getter=lambda: M1.tau, units='hr', element=M1),
        bst.Metric(name='tau2', getter=lambda: M2.tau, units='hr', element=M2),
        ]
    
    np.random.seed(3221)
    samples = model.sample(100, rule='L')
    model.load_samples(samples)
    model.evaluate()
    
    D, p = model.kolmogorov_smirnov_d(thresholds=[1, 1.5]) # Just make sure it works for now
    # TODO: Add tests that make sense for comparing statistics
    
if __name__ == '__main__':
    test_parameter_hook()
    test_pearson_r()
    test_spearman_r()
    test_kendall_tau()
    test_model_index()
    test_model_sample()
    test_copy()
    test_model_exception_hook()
    test_kolmogorov_smirnov_d()