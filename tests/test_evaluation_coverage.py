# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import biosteam as bst
import numpy as np
import pandas as pd
import os
import tempfile
from numpy.testing import assert_allclose
from chaospy import distributions as shape

def create_model(ID):
    bst.F.set_flowsheet(ID)
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=1)
    with bst.System(ID) as system:
        M1 = bst.Mixer('M1', feed, 'product')
    model = bst.Model(system)

    @model.parameter(element=feed, units='kmol/hr', bounds=(5, 15),
                     distribution='uniform', coupled=True)
    def set_water(water):
        feed.imol['Water'] = water

    @model.parameter(element=feed, units='kmol/hr', bounds=(0.5, 2),
                     baseline=1, distribution='triangular')
    def set_ethanol(ethanol):
        feed.imol['Ethanol'] = ethanol

    @model.indicator(units='kg/hr', element=M1)
    def product_flow():
        return M1.outs[0].F_mass

    @model.indicator
    def ethanol_fraction():
        return M1.outs[0].imass['Ethanol'] / M1.outs[0].F_mass

    return model

def test_model_properties_and_errors():
    model = create_model('model_properties_and_errors')
    assert len(model) == 2
    assert model.specification is None
    model.specification = lambda: model.system.simulate()
    assert callable(model.specification)
    model.specification = None
    with pytest.raises(AttributeError):
        model.specification = 1
    baseline = model.get_baseline_scenario(array=True)
    assert_allclose(baseline, [10, 1], rtol=1e-2)
    parameters = model.parameters
    model.parameters = parameters
    assert model.parameters == parameters
    model.optimized_parameters = ()
    assert model.optimized_parameters == ()
    model.metrics = model.indicators
    assert len(model._metrics) == 2
    model._metrics = model._metrics
    indicators = model.indicators
    with pytest.raises(ValueError):
        model.indicators = [1]
    model.indicators = indicators
    with pytest.raises(ValueError):
        model.sample(3, 'invalid')
    with pytest.raises(ValueError):
        model.exception_hook = 'invalid'
    with pytest.raises(ValueError):
        model.exception_hook = 1
    model.exception_hook = 'ignore'
    assert model.exception_hook(None, None) is None
    model.exception_hook = 'raise'
    with pytest.raises(RuntimeError):
        model.exception_hook(RuntimeError(), [1, 2])
    with pytest.raises(TypeError):
        model.load_samples([1, 2])
    with pytest.raises(ValueError):
        model.load_samples(np.ones([2, 2, 2]))
    with pytest.raises(ValueError):
        model.load_samples(np.ones([2, 3]))
    with pytest.raises(RuntimeError):
        model.evaluate()
    with pytest.raises(ValueError):
        model.parameter('flow')
    p = model.parameter(lambda x: None, name='No baseline')
    with pytest.raises(RuntimeError):
        model.get_baseline_scenario()
    model.parameters = model.parameters[:-1]
    with pytest.raises(ValueError):
        model._load_sample_order(np.ones([3, 2]), model.parameters, None, 'invalid')
    assert 'parameters' in repr(model)
    assert 'Model' in model._repr(None)
    model.show(p=1, m=1)
    empty_model = bst.Model(None)
    empty_model.show()
    assert empty_model.get_distribution_summary() is None

def test_model_parameter_and_indicator_redefinition():
    model = create_model('model_parameter_and_indicator_redefinition')
    set_water, set_ethanol = model.parameters
    product_flow, ethanol_fraction = model.indicators
    p = model.parameter(set_water) # Replaces old parameter
    assert p.name == set_water.name
    assert len(model.parameters) == 2
    with pytest.raises(ValueError):
        model.parameter(set_water, safe=True)
    p = model.parameter(lambda x: None, element=set_water.element, name='Other')
    assert len(model.parameters) == 3
    p = model.optimized_parameter(lambda x: None, element=set_water.element, name='Optimized', bounds=(0, 1))
    assert len(model.optimized_parameters) == 1
    p = model.optimized_parameter(lambda x: None, element=set_water.element, name='Optimized', bounds=(0, 1))
    assert len(model.optimized_parameters) == 1
    indicator = model.indicator(product_flow)
    assert indicator.name == product_flow.name
    assert len(model.indicators) == 2
    model_copy = model.copy()
    assert len(model_copy.parameters) == len(model.parameters)

def test_mock_variables():
    model = create_model('mock_variables')
    set_water = model.parameters[0]
    product_flow = model.indicators[0]
    mock_parameter = set_water.mockup()
    assert 'MockVariable' in repr(mock_parameter)
    p = model.parameter(lambda x: None, element=mock_parameter.element, name='Other water')
    mock_parameter = bst.MockVariable('Other water 2', 'kmol/hr', set_water.element)
    p = model.parameter(mock_parameter, element=None)
    assert p is not None
    mock_indicator = bst.MockVariable('Other flow', 'kg/hr', product_flow.element)
    indicator = model.indicator(mock_indicator)
    assert indicator.name == 'Other flow'
    assert indicator.units == 'kg/hr'

def test_variable_descriptions():
    model = create_model('variable_descriptions')
    set_water, set_ethanol = model.parameters
    product_flow, ethanol_fraction = model.indicators
    assert set_water.name_with_units == 'Water [kmol/hr]'
    assert ethanol_fraction.name_with_units == 'Ethanol fraction'
    assert isinstance(set_water.short_description, str)
    assert isinstance(set_water.short_description_with_units, str)
    long = bst.Variable('A very long variable name used for testing descriptions', 'kg', None)
    assert len(long.short_description) <= 31
    assert len(long.short_description_with_units) <= 31
    assert 'kg' in long.label()
    assert 'Water' in set_water.describe(format_units=True)
    assert '[' in set_water.describe(distribution=False, bounds=False)
    assert '(' in set_ethanol.describe(distribution=False)
    nounits = bst.Parameter('value', lambda x: None, None, None, None,
                            None, 0.5, (0, 1), False, None, None)
    assert nounits.describe() == 'Value (Uniform; 0, 1)'
    assert nounits.describe(distribution=False) == 'Value (0, 0.5, 1)'
    assert 'Value' in repr(nounits)
    nounits.show()
    with pytest.raises(ValueError):
        bst.Variable.check_indices_unique([set_water, set_water])

def test_parameter():
    model = create_model('parameter')
    set_water, set_ethanol = model.parameters
    feed = model.system.flowsheet.feed
    M1 = model.system.flowsheet.M1
    assert set_water.unit is M1
    assert set_water.subsystem is not None
    set_water(12)
    assert_allclose(M1.outs[0].imol['Water'], 12, rtol=1e-2)
    set_ethanol.hook = lambda x: 2 * x
    set_ethanol(0.5)
    assert_allclose(feed.imol['Ethanol'], 1, rtol=1e-2)
    set_ethanol.hook = None
    p = bst.Parameter('M1 parameter', lambda x: None, M1, None, None,
                      None, None, (0, 1), True, None, None)
    assert p.subsystem is None
    with pytest.raises(RuntimeError):
        p.simulate()
    p = bst.Parameter('M1 parameter', lambda x: None, None, None,
                      shape.Uniform(0, 1), None, None, None, False, None, None)
    assert p.bounds == (0, 1)
    p.simulate()
    with pytest.raises(ValueError):
        bst.Parameter('Bad', lambda x: None, None, None, 'normal', None,
                      None, (0, 1), False, None, None)
    assert p._info().startswith('Parameter')
    parameters = list(model.parameters)
    bst.Parameter.sort_parameters(parameters)
    assert len(parameters) == 2
    bst.Parameter.sort_parameters([])
    bst.Parameter.sort_parameters([p])
    with pytest.raises(ValueError):
        bst.Parameter.sort_parameters([p, set_water])

def test_parameters_from_df_errors():
    model = create_model('parameters_from_df_errors')
    df = pd.DataFrame({
        'Parameter name': ['Water flow'],
        'Element': ['feed'],
        'Kind': ['coupled'],
        'Units': ['kmol/hr'],
        'Baseline': [10],
        'Shape': ['Uniform'],
        'Lower': [5],
        'Midpoint': [10],
        'Upper': [15],
        'Load statement': ['feed.imol["Water"] = x'],
    })
    with pytest.raises(ValueError):
        model.parameters_from_df(df)
    df['Midpoint'] = [float('nan')]
    with tempfile.TemporaryDirectory() as tmp:
        file = os.path.join(tmp, 'parameters.csv')
        df.to_csv(file, index=False)
        model.parameters_from_df(file)
    assert len(model.parameters) == 3

def test_distribution_summary():
    model = create_model('distribution_summary')
    with tempfile.TemporaryDirectory() as tmp:
        file = os.path.join(tmp, 'summary.xlsx')
        tables = model.get_distribution_summary(file)
        assert os.path.exists(file)
    assert set(tables) == {'Uniform', 'Triangle'}

def test_model_evaluate_options():
    model = create_model('model_evaluate_options')
    samples = model.sample(4, 'L', seed=0)
    with tempfile.TemporaryDirectory() as tmp:
        samples_file = os.path.join(tmp, 'folder', 'samples.pkl')
        model.load_samples(samples, file=samples_file, autosave=True)
        model.load_samples(samples, file=samples_file, autoload=True)
        file = os.path.join(tmp, 'other_folder', 'results.pkl')
        model.evaluate(notify=2, file=file, autosave=1)
        table = model.table.copy()
        model.evaluate(file=file, autoload=True, notify=1)
        assert_allclose(model.table.values, table.values, rtol=1e-2)
        model.load_pickled_results(file)
        assert_allclose(model.table.values, table.values, rtol=1e-2)
        model.load_samples(samples[:1])
        model.evaluate(file=file, autoload=True)
    assert np.isfinite(model.table.values).all()
    values = model(model.get_baseline_scenario())
    assert values.shape == (2,)
    set_water = model.parameters[0]
    set_water.active = False
    model.load_samples(samples)
    model.evaluate()
    set_water.active = True
    assert np.isfinite(model.table.values).all()

def test_model_exception_hook_and_retry():
    model = create_model('model_exception_hook_and_retry')
    state = {'fail': True}
    def indicator():
        if state['fail']:
            state['fail'] = False
            raise RuntimeError('failed')
        return 1.
    model.indicator(indicator, name='Failing')
    model.retry_evaluation = False
    model.exception_hook = lambda exception, sample: [1., 1., 1.]
    model.load_samples(model.sample(2, 'L', seed=0))
    model.evaluate()
    assert np.isfinite(model.table.values).all()
    state['fail'] = True
    model.exception_hook = lambda exception, sample: 1.
    with pytest.raises(RuntimeError):
        model.evaluate()
    def specification():
        if state['fail']:
            state['fail'] = False
            raise RuntimeError('failed')
        model.system.simulate()
    state['fail'] = True
    model.specification = specification
    model.retry_evaluation = True
    model.exception_hook = None
    model.evaluate()
    assert np.isfinite(model.table.values).all()

def test_single_point_sensitivity():
    model = create_model('single_point_sensitivity')
    model.parameters[1].hook = lambda x: x
    baseline, lb, ub = model.single_point_sensitivity()
    assert lb.shape == ub.shape == (2, 2)
    assert np.isfinite(lb.values).all() and np.isfinite(ub.values).all()
    baseline, lb, ub = model.single_point_sensitivity(array=True)
    assert lb.shape == (2, 2)
    count = [0]
    def evaluate(sample, **kwargs):
        count[0] += 1
        return [count[0]]
    baseline, lb, ub = model.single_point_sensitivity(
        evaluate=evaluate, indicators=model.indicators[:1], array=True,
    ) # Inconsistent model is printed
    assert lb.shape == (2, 1)

def test_model_optimize():
    model = create_model('model_optimize')
    set_water, set_ethanol = model.parameters
    product_flow, ethanol_fraction = model.indicators
    result = model.optimize(lambda: -ethanol_fraction(), parameters=[set_ethanol],
                            method='SLSQP')
    assert_allclose(result.x, [2], rtol=5e-2)
    result = model.optimize(lambda: -ethanol_fraction(), parameters=[set_water],
                            method='differential evolution',
                            optimizer_options=dict(maxiter=2, popsize=2, seed=0),
                            convergence_model='average')
    assert result.convergence_model is not None
    with pytest.raises(ValueError):
        model.optimize(lambda: 0, method='invalid')

def test_evaluate_across_coordinate():
    model = create_model('evaluate_across_coordinate')
    set_water, set_ethanol = model.parameters
    feed = model.system.flowsheet.feed
    model.load_samples(model.sample(3, 'L', seed=0))
    with tempfile.TemporaryDirectory() as tmp:
        file = os.path.join(tmp, 'coordinate.xlsx')
        data = model.evaluate_across_coordinate(
            'Water', set_water, [5, 10], xlfile=file,
        )
        assert os.path.exists(file)
        file = os.path.join(tmp, 'multi_coordinate.xlsx')
        data = model.evaluate_across_coordinate(
            ('Water', 'Ethanol'), lambda x, y: None, [(5, 1), (10, 2)],
            xlfile=file, multi_coordinate=True, notify_coordinate=False,
        )
    assert all([i.shape == (3, 2) for i in data.values()])
    def set_price(x): feed.price = x
    model.indicator(lambda: feed.price, name='Price')
    data = model.evaluate_across_coordinate(
        'Price', set_price, [1, 2], simulation_independent_coordinate=True,
    )
    assert_allclose(data[model.indicators[-1].index], [[1, 2]] * 3, rtol=1e-2)
    with pytest.raises(ValueError):
        model.evaluate_across_coordinate(
            'Price', set_price, [1, 2], simulation_independent_coordinate=True,
            f_evaluate=model.evaluate,
        )
    model._samples = None
    with pytest.raises(RuntimeError):
        model.evaluate_across_coordinate(
            'Price', set_price, [1, 2], simulation_independent_coordinate=True,
        )

def test_model_correlations():
    model = create_model('model_correlations')
    model.load_samples(model.sample(5, 'L', seed=0))
    model.evaluate()
    with pytest.warns(DeprecationWarning):
        rho = model.spearman()
    assert rho.shape == (2, 2)
    for filter in ('omit nan', 'raise nan', 'none', lambda x, y: (x, y)):
        rho, p = model.pearson_r(filter=filter)
        assert rho.shape == (2, 2)
    with pytest.raises(ValueError):
        model.kolmogorov_smirnov_d(thresholds=[1.])
    model.table.iloc[0, -1] = np.nan
    rho, p = model.pearson_r(filter='propagate nan')
    assert np.isnan(rho.values[:, -1]).all()
    with pytest.raises(ValueError):
        model.pearson_r(filter='raise nan')

def test_convergence_model_options():
    bst.F.set_flowsheet('convergence_model_options')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=10)
    recycle = bst.Stream('recycle')
    with bst.System('convergence_model_options') as system:
        M1 = bst.Mixer('M1', [feed, recycle])
        S1 = bst.Splitter('S1', M1-0, ['product', recycle], split=0.5)
    model = bst.Model(system)

    @model.parameter(bounds=(5, 15), coupled=True)
    def set_water(water):
        feed.imol['Water'] = water

    @model.parameter(bounds=(0.3, 0.6), coupled=True)
    def set_split(split):
        S1.split[:] = split

    @model.parameter(bounds=(0, 1))
    def set_price(price):
        feed.price = price

    @model.indicator
    def product_flow():
        return S1.outs[0].F_mol

    model.load_samples(model.sample(4, 'L', seed=0), sort=True)
    parameters = model.parameters
    for kwargs in (
            dict(model_type='linear regressor', save_prediction=True),
            dict(model_type='intercept linear regressor', interaction_pairs=True),
            dict(model_type='average', save_prediction=True, local_weighted=False, recess=1),
            dict(model_type=bst.InterceptLinearRegressor, local_weighted=False, nfits=1),
            dict(model_type=bst.LinearRegressor, local_weighted=False, save_prediction=True, nfits=1),
        ):
        convergence_model = bst.ConvergenceModel(parameters, **kwargs)
        model.evaluate(convergence_model=convergence_model)
        assert np.isfinite(model.table.values).all()
        assert len(convergence_model.responses) > 0
    convergence_model.R2_predicted(last=2)
    repr(bst.LinearRegressor()), repr(bst.InterceptLinearRegressor()), repr(bst.Average())
    for response in convergence_model.responses:
        assert str(response)
        assert repr(response)
    model.evaluate(convergence_model='average')
    assert np.isfinite(model.table.values).all()
    with pytest.raises(ValueError):
        bst.ConvergenceModel(parameters, model_type='invalid')
    with pytest.raises(ValueError):
        bst.ConvergenceModel(parameters, local_weighted=True, recess=1)
    with pytest.raises(ValueError):
        bst.ConvergenceModel(parameters, local_weighted=True, nfits=1)
    with pytest.raises(ValueError):
        bst.ConvergenceModel([])
    with pytest.raises(ValueError):
        bst.NullConvergenceModel([])
    other = bst.Parameter('other', lambda x: None, None, bst.System(None, ()), None,
                          None, None, (0, 1), True, None, None)
    with pytest.raises(ValueError):
        bst.ConvergenceModel([parameters[0], other])
    with pytest.raises(ValueError):
        bst.NullConvergenceModel([parameters[0], other])
    convergence_model = bst.NullConvergenceModel(parameters)
    model.evaluate(convergence_model=convergence_model)
    assert np.isfinite(model.table.values).all()

def test_generic_response():
    bst.F.set_flowsheet('generic_response')
    bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    feed = bst.Stream('feed', Water=10, Ethanol=10)
    response = bst.GenericResponse(feed, 'T', bst.Average(), 'K')
    response.update_limits(300)
    response.update_limits(400)
    response.set(500)
    assert_allclose(feed.T, 400, rtol=1e-2)
    response.set(200)
    assert_allclose(response.get(), 300, rtol=1e-2)
    response.set(350)
    assert_allclose(response.get(), 350, rtol=1e-2)
    assert 'GenericResponse' in repr(response)
    assert hash(response) == hash((feed, 'T'))
    response.update_limits(450)
    response.update_limits(250)
    assert response.max == 450 and response.min == 250
    response.parameters = [0]
    X = np.array([[1.], [2.], [3.]])
    y = np.array([300., 310., 320.])
    response.fit(X, y)
    assert_allclose(response.predict(np.array([2.])), 310, rtol=1e-2)
    value = response.predict_locally(np.array([2.]), X, y, 'cityblock', lambda x: np.exp(-x))
    assert_allclose(value, 310, rtol=1e-2)
    value = response.predict_locally(np.array([2.5]), X, y, 'cityblock', lambda x: np.exp(-x))
    assert 250 <= value <= 450

if __name__ == '__main__':
    test_model_properties_and_errors()
    test_model_parameter_and_indicator_redefinition()
    test_mock_variables()
    test_variable_descriptions()
    test_parameter()
    test_parameters_from_df_errors()
    test_distribution_summary()
    test_model_evaluate_options()
    test_model_exception_hook_and_retry()
    test_single_point_sensitivity()
    test_model_optimize()
    test_evaluate_across_coordinate()
    test_model_correlations()
    test_convergence_model_options()
    test_generic_response()
