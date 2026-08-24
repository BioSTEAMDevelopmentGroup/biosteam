# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Copyright (C) 2026-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Tests that the biosteam <-> hensmith circular dependency is import-safe in
both directions and that biosteam lazily re-exports HeatExchangerNetwork
(and, with a DeprecationWarning, the network synthesis helpers) from
hensmith via PEP 562, so that `import biosteam` never requires hensmith to
be initialized first and vice versa.
"""
import subprocess
import sys
import pytest

def _run(code):
    # Import-ordering semantics can only be tested in a cold interpreter;
    # capture output so a child failure surfaces its traceback instead of
    # an opaque "returned non-zero exit status 1".
    result = subprocess.run(
        [sys.executable, '-c', code], capture_output=True, text=True,
    )
    if result.returncode:
        raise AssertionError(
            f"subprocess failed with exit code {result.returncode}\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )

def test_biosteam_first_import_order():
    _run(
        "import biosteam as bst\n"
        "import hensmith\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        # Supported lazy names are cached on first access (PEP 562 shim
        # memoization) and reported by the companion __dir__.
        "assert 'HeatExchangerNetwork' in vars(bst)\n"
        "assert 'HeatExchangerNetwork' in dir(bst)\n"
        "assert 'HeatExchangerNetwork' in dir(bst.facilities)\n"
        # The from-import and star-import forms (folded from a separate
        # subprocess; they add no ordering coverage of their own).
        "from biosteam import HeatExchangerNetwork\n"
        "from biosteam.facilities import HeatExchangerNetwork as HXN2\n"
        "assert HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert HXN2 is hensmith.HeatExchangerNetwork\n"
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
        # Star-importing the facilities subpackage must NOT provide the name:
        # keeping it out of facilities.__all__ is the guard that stops
        # biosteam/__init__'s star-import of facilities from re-creating the
        # biosteam <-> hensmith import cycle.
        "ns = {}\n"
        "exec('from biosteam.facilities import *', ns)\n"
        "assert 'HeatExchangerNetwork' not in ns\n"
    )

def test_hensmith_first_import_order():
    _run(
        "import hensmith\n"
        "import biosteam as bst\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        # `from biosteam import *` resolves HeatExchangerNetwork eagerly
        # (it is in biosteam.__all__), which is safe only while no module
        # executed during hensmith's initialization star-imports biosteam
        # at module scope — exercise it under hensmith-first ordering too.
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
    )

def test_star_import_initializes_hensmith_from_scratch():
    # `from biosteam import *` in a fresh interpreter, hensmith never
    # imported: the star-import itself must initialize hensmith through the
    # lazy re-export — hensmith's own `import biosteam` then finds the
    # already-completed module. This is the common real-world consumer path.
    _run(
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "import hensmith\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
    )

def test_missing_hensmith_degrades_to_AttributeError():
    # With hensmith unavailable, biosteam must still import, and attribute
    # access must fail with AttributeError — not leak ModuleNotFoundError
    # through hasattr/getattr — with a message pointing at hensmith.
    _run(
        "import sys\n"
        "sys.modules['hensmith'] = None\n"  # makes `import hensmith` raise ImportError
        "import biosteam as bst\n"
        "assert not hasattr(bst, 'HeatExchangerNetwork')\n"
        "assert not hasattr(bst.facilities, 'HeatExchangerNetwork')\n"
        "assert getattr(bst, 'HeatExchangerNetwork', None) is None\n"
        "for module in (bst, bst.facilities):\n"
        "    try:\n"
        "        module.HeatExchangerNetwork\n"
        "    except AttributeError as e:\n"
        "        assert 'hensmith' in str(e)\n"
        "    else:\n"
        "        raise AssertionError('AttributeError not raised')\n"
    )

def test_HeatExchangerNetwork_not_in_facilities_all():
    # Negative guard: re-adding the name to facilities.__all__ would make
    # biosteam/__init__'s star-import of facilities resolve it eagerly and
    # re-create the biosteam <-> hensmith import cycle. Fail loudly if a
    # future "fix" tries.
    import biosteam as bst
    assert 'HeatExchangerNetwork' not in bst.facilities.__all__
    for name in bst.facilities._HENSMITH_DEPRECATED_NAMES:
        assert name not in bst.facilities.__all__
        assert name not in bst.__all__

def test_deprecated_helpers_forward_with_warning():
    import biosteam as bst
    import hensmith
    for module in (bst, bst.facilities):
        for name in bst.facilities._HENSMITH_DEPRECATED_NAMES:
            with pytest.warns(DeprecationWarning, match='hensmith'):
                value = getattr(module, name)
            assert value is getattr(hensmith, name)
            # Deprecated aliases are not cached, so every access warns.
            assert name not in vars(module)

def test_create_all_facilities_constructs_hensmith_network():
    # The real downstream seam: create_facilities instantiates
    # bst.HeatExchangerNetwork (HXN=True by default) through the lazy
    # re-export inside a MockSystem context.
    import biosteam as bst
    import hensmith
    bst.settings.set_thermo(['Water'], cache=True)
    bst.main_flowsheet.set_flowsheet('hensmith_integration_HXN_smoke')
    try:
        units = bst.create_all_facilities(
            CT=False, CWP=False, CIP=False, FWT=False, ADP=False,
            WWT=False, CHP=False, PWC=False,  # HXN=True by default
        )
        hxns = [i for i in units if isinstance(i, hensmith.HeatExchangerNetwork)]
        assert len(hxns) == 1
    finally:
        bst.main_flowsheet.clear()

if __name__ == '__main__':
    test_biosteam_first_import_order()
    test_hensmith_first_import_order()
    test_star_import_initializes_hensmith_from_scratch()
    test_missing_hensmith_degrades_to_AttributeError()
    test_HeatExchangerNetwork_not_in_facilities_all()
    test_deprecated_helpers_forward_with_warning()
    test_create_all_facilities_constructs_hensmith_network()
