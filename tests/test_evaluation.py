# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest

def test_model():
    from biosteam.examples import ethanol_subsystem as ethanol_sys
    from biosteam import main_flowsheet as F
    import biosteam as bst
    
    # Make sure repeated metrics raise an error
    biorefinery = bst.process_tools.UnitGroup('Biorefinery', ethanol_sys.units)
    heating_duty = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    heating_duty_repeated = bst.Metric('Heating duty', biorefinery.get_heating_duty, 'GJ/hr')
    with pytest.raises(ValueError):
        model = bst.Model(ethanol_sys, [heating_duty, heating_duty_repeated])
    
    model = bst.Model(ethanol_sys, [heating_duty])
    
    # Make sure repeated parameters raise an error
    R301 = F.flowsheet.ethanol_subsystem_example.unit.R301
    @model.parameter(element=R301)
    def set_efficiency(efficiency):
        R301.efficiency = efficiency
    
    with pytest.raises(ValueError):
        @model.parameter(element=R301)
        def set_efficiency(efficiency):
            R301.efficiency = efficiency
    
    bst.process_tools.default()
    
def test_model_exception_hook():
    import biosteam as bst
    import pytest
    from biorefineries import lipidcane as lc
    from chaospy import distributions as shape
    from warnings import simplefilter
    import numpy as np
    bst.settings.set_thermo(lc.chemicals)
    simplefilter("ignore")
    IRR_metric = bst.Metric('Internal rate of return', lc.lipidcane_tea.solve_IRR)
    metrics = [IRR_metric]
    lipidcane_model = bst.Model(lc.lipidcane_sys, metrics)
    baseline = lc.lipidcane.F_mass
    distribution = shape.Triangle(-baseline , baseline , 2*baseline) # Negative value should fail
    
    @lipidcane_model.parameter(element=lc.lipidcane, distribution=distribution, units='kg/hr')
    def set_lipidcane_flow_rate(flow_rate):
        lc.lipidcane.F_mass = flow_rate
    
    samples = lipidcane_model.sample(15, 'L')
    lipidcane_model.load_samples(samples)
    
    # Without an exception hook, the same behavior will result (NaN values for failed evaluations)
    lipidcane_model.evaluate()
    assert np.isnan(lipidcane_model.table.values).any()
    
    InfeasibleRegion = bst.exceptions.InfeasibleRegion
    
    # This will provide a more understandable IRR result for infeasible regions
    def exception_hook(exception, sample): 
        if isinstance(exception, InfeasibleRegion):
            return [0]
        else: raise exception
    lipidcane_model.exception_hook = exception_hook
    lipidcane_model.evaluate()
    assert not np.isnan(lipidcane_model.table.values).any()
    
    # This will raise an exception due to negative flow rates
    def exception_hook(exception, sample): 
        if isinstance(exception, InfeasibleRegion):
            raise exception
    lipidcane_model.exception_hook = exception_hook
    with pytest.raises(InfeasibleRegion): lipidcane_model.evaluate()
    
    # This will raise an exception regardless
    def exception_hook(exception, sample): 
        raise exception
    lipidcane_model.exception_hook = exception_hook
    with pytest.raises(InfeasibleRegion): lipidcane_model.evaluate()
    
    # Here is another cool thing we could do in the case where 
    # some metrics are expected to fail
    bad_metric = bst.Metric('bad metric', lambda: 1/0)
    lipidcane_model.metrics = (IRR_metric, bad_metric)
    lipidcane_model.load_samples(samples) # Metrics changed, so need to reload sample
    def exception_hook(exception, sample): 
        if not isinstance(exception, ZeroDivisionError): return
        lc.lipidcane_sys.simulate()
        values = []
        for i in lipidcane_model.metrics:
            try: x = i()
            except: x = None
            values.append(x)
        return values
    lipidcane_model.exception_hook = exception_hook
    lipidcane_model.evaluate()
    bad_metric_results = lipidcane_model.table[bad_metric.index]
    IRR_metric_results = lipidcane_model.table[IRR_metric.index]
    assert np.isnan(bad_metric_results).all()
    assert not np.isnan(IRR_metric_results).all()
    
    # Note that the possibilities are infinite here...
    