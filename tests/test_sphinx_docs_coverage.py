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
import os
import io
import tempfile
import importlib
sd = importlib.import_module("biosteam.utils.sphinx_docs")

def example_function(x):
    """Return x."""
    return x

def test_autodoc_directive():
    text = sd.autodoc_directive(bst.Mixer, inherited_members=True)
    assert text.startswith('.. autoclass:: ')
    assert ':inherited-members:' in text
    text = sd.autodoc_directive(RuntimeError, members=False)
    assert text.startswith('.. autoexception:: ')
    assert ':members:' not in text
    text = sd.autodoc_directive(sd)
    assert text.startswith('.. automodule:: ')
    text = sd.autodoc_directive(example_function)
    assert text.startswith('.. autofunction:: ')
    text = sd.autodoc_directive(bst.Mixer.simulate)
    assert text.startswith('.. automethod:: ')
    with pytest.raises(TypeError):
        sd.autodoc_directive(1)

def test_sphinx_docs():
    pytest.importorskip('sphinx')
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.raises(FileNotFoundError):
            sd.sphinx_docs('text', confdir=tmp)
    if not os.path.isfile(os.path.join(sd.default_confdir, 'conf.py')):
        pytest.skip('docs folder not available')
    with pytest.raises(ValueError):
        sd.sphinx_docs('text', format='txt')
    text = sd.sphinx_docs(
        example_function, warnings=io.StringIO(),
        confoverrides={'extensions': ['sphinx.ext.autodoc']},
    )
    assert 'example_function' in text

if __name__ == '__main__':
    test_autodoc_directive()
    test_sphinx_docs()
