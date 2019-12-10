# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:52:11 2019

@author: yoelr
"""

__all__ = ('assert_unit_results', 'assert_unit_streams')

def assert_unit_results(unit, unit_results):
    assert str(unit.results()) == unit_results, "significant changes in simulation results"

def assert_stream_results(unit, stream_results):
    assert unit._info('K', 'Pa', 'kmol/hr', False, 1000) == stream_results, "significant changes in stream results"