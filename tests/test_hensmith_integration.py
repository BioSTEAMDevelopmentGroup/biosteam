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
both directions: biosteam imports hensmith at the very end of its own
initialization, and hensmith binds HeatExchangerNetwork into biosteam and
biosteam.facilities when it finishes initializing, so the name is available
eagerly whichever package is imported first.
"""
import subprocess
import sys

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
        "import sys\n"
        "import biosteam as bst\n"
        # Eager: importing biosteam initializes hensmith and binds the name
        # without any further import on the user's part.
        "assert 'hensmith' in sys.modules\n"
        "assert 'HeatExchangerNetwork' in vars(bst)\n"
        "assert 'HeatExchangerNetwork' in vars(bst.facilities)\n"
        "import hensmith\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        # The from-import and star-import forms.
        "from biosteam import HeatExchangerNetwork\n"
        "from biosteam.facilities import HeatExchangerNetwork as HXN2\n"
        "assert HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert HXN2 is hensmith.HeatExchangerNetwork\n"
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
        # Star-importing the facilities subpackage must NOT provide the name:
        # biosteam/__init__ star-imports facilities before hensmith can bind
        # it, so listing it in facilities.__all__ would break `import biosteam`.
        "ns = {}\n"
        "exec('from biosteam.facilities import *', ns)\n"
        "assert 'HeatExchangerNetwork' not in ns\n"
    )

def test_hensmith_first_import_order():
    # biosteam's `import hensmith` runs while hensmith is still initializing
    # (HeatExchangerNetwork not yet defined); hensmith must bind the name
    # itself once it finishes.
    _run(
        "import hensmith\n"
        "import biosteam as bst\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
    )

def test_star_import_initializes_hensmith_from_scratch():
    # `from biosteam import *` in a fresh interpreter, hensmith never
    # imported by the user: the common real-world consumer path.
    _run(
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "import hensmith\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
    )

def test_missing_hensmith_degrades_gracefully():
    # With hensmith unavailable, biosteam must still import (with the name
    # absent from the namespace and from __all__, so star-imports keep
    # working) rather than fail at import time.
    _run(
        "import sys\n"
        "sys.modules['hensmith'] = None\n"  # makes `import hensmith` raise ModuleNotFoundError
        "import biosteam as bst\n"
        "assert not hasattr(bst, 'HeatExchangerNetwork')\n"
        "assert not hasattr(bst.facilities, 'HeatExchangerNetwork')\n"
        "assert 'HeatExchangerNetwork' not in bst.__all__\n"
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "assert 'HeatExchangerNetwork' not in ns\n"
    )

def test_broken_hensmith_is_not_swallowed():
    # Only a missing hensmith is tolerated; an error raised while hensmith
    # itself initializes must propagate out of `import biosteam`.
    _run(
        "import sys, types\n"
        "class BrokenFinder:\n"
        "    @staticmethod\n"
        "    def find_spec(name, path=None, target=None):\n"
        "        if name == 'hensmith':\n"
        "            import importlib.util\n"
        "            loader = types.SimpleNamespace(\n"
        "                create_module=lambda spec: None,\n"
        "                exec_module=lambda module: exec('raise RuntimeError(\"hensmith is broken\")'),\n"
        "            )\n"
        "            return importlib.util.spec_from_loader(name, loader)\n"
        "sys.meta_path.insert(0, BrokenFinder)\n"
        "try:\n"
        "    import biosteam\n"
        "except RuntimeError as e:\n"
        "    assert 'hensmith is broken' in str(e)\n"
        "else:\n"
        "    raise AssertionError('error inside hensmith was swallowed')\n"
    )

def test_HeatExchangerNetwork_not_in_facilities_all():
    # Negative guard: biosteam/__init__ star-imports facilities before
    # hensmith binds the name, so re-adding it to facilities.__all__ would
    # break `import biosteam`. Fail loudly if a future "fix" tries.
    import biosteam as bst
    assert 'HeatExchangerNetwork' not in bst.facilities.__all__
    assert 'HeatExchangerNetwork' in bst.__all__

def test_create_all_facilities_constructs_hensmith_network():
    # The real downstream seam: create_facilities instantiates
    # bst.HeatExchangerNetwork (HXN=True by default) inside a MockSystem
    # context.
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
    test_missing_hensmith_degrades_gracefully()
    test_broken_hensmith_is_not_swallowed()
    test_HeatExchangerNetwork_not_in_facilities_all()
    test_create_all_facilities_constructs_hensmith_network()
