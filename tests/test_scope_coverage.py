# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2024, Joy Zhang <joycheung1994@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import pytest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import biosteam as bst
import numpy as np
import pandas as pd
import os
import tempfile
from types import SimpleNamespace
from numpy.testing import assert_allclose
from biosteam.utils import Scope, SystemScope

def create_subject(ID):
    subject = SimpleNamespace(ID=ID, state=np.array([1., 2.]), T=np.array([300.]))
    subject.scope = Scope(subject, ['state', 'T'])
    return subject

def record(subject, scope, times):
    for t in times:
        subject.state = subject.state + 1
        subject.T = subject.T + 10
        scope(t)

def test_scope():
    subject = create_subject('U1')
    scope = subject.scope
    assert repr(scope) == '<Scope: U1>'
    record(subject, scope, [0., 1., 2.])
    assert scope.record.shape == (3, 3)
    assert_allclose(scope.time_series, [0, 1, 2])
    scope.header = None
    assert scope.header == [('U1', 'state_0'), ('U1', 'state_1'), ('U1', 'T_0')]
    header = [('U1', 'a'), ('U1', 'b'), ('U1', 'c')]
    scope.header = header
    assert scope.header == header
    with pytest.raises(ValueError):
        scope.header = [('U1', 'a', 'x'), ('U1', 'b', 'x'), ('U1', 'c', 'x')]
    fig, ax = scope.plot_time_series('T')
    fig, ax = scope.plot_time_series('state')
    plt.close('all')
    scope.pop()
    assert scope.record.shape == (2, 3)
    scope.reset_cache()
    assert len(scope.time_series) == 0

def test_system_scope():
    bst.F.set_flowsheet('system_scope')
    system = bst.System('system_scope', ())
    system._state_header = None
    U1 = create_subject('U1')
    U2 = create_subject('U2')
    U2.scope.header = None
    U1.scope.header = pd.MultiIndex.from_tuples([('U1', 'a'), ('U1', 'b'), ('U1', 'c')])
    scope = SystemScope(system, U1, U2)
    assert repr(scope) == '<SystemScope: system_scope>'
    assert list(scope.sol_header) == [('-', 't [d]')]
    for t in (0., 1., 1., 2., 3.): # Repeated time is replaced
        for subject in (U1, U2):
            subject.state = subject.state + 1
            subject.T = subject.T + 10
        scope(t)
    assert_allclose(scope.time_series, [0, 1, 2, 3])
    df = scope.export()
    assert df.shape == (4, 7)
    df = scope.export(t_eval=[0.5, 1.5])
    assert df.shape == (2, 7)
    assert_allclose(df.iloc[:, 0], [0.5, 1.5])
    with pytest.raises(RuntimeError):
        scope.export(t_eval=[10.])
    with tempfile.TemporaryDirectory() as tmp:
        for ext in ('npy', 'csv', 'tsv', 'xlsx'):
            file = os.path.join(tmp, 'data.' + ext)
            scope.export(file)
            assert os.path.exists(file)
        with pytest.raises(ValueError):
            scope.export(os.path.join(tmp, 'data.txt'))
    scope.sol_header = [('-', 't [d]')]
    assert scope.sol_header == [('-', 't [d]')]
    scope.reset_cache()
    assert len(scope.time_series) == 0
    assert scope.sol is None
    with pytest.raises(AttributeError):
        SystemScope(system, SimpleNamespace())
    with pytest.raises(TypeError):
        SystemScope(system, SimpleNamespace(scope=1))

if __name__ == '__main__':
    test_scope()
    test_system_scope()
