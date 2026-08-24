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
from hensmith (PEP 562), so that `import biosteam` never requires hensmith
to be initialized first and vice versa.
"""
import subprocess
import sys

def _run(code):
    subprocess.run([sys.executable, '-c', code], check=True)

def test_biosteam_first_import_order():
    _run(
        "import biosteam as bst\n"
        "import hensmith\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
    )

def test_hensmith_first_import_order():
    _run(
        "import hensmith\n"
        "import biosteam as bst\n"
        "assert bst.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert bst.facilities.HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
    )

def test_from_imports_and_star_import():
    _run(
        "from biosteam import HeatExchangerNetwork\n"
        "from biosteam.facilities import HeatExchangerNetwork as HXN2\n"
        "ns = {}\n"
        "exec('from biosteam import *', ns)\n"
        "import hensmith\n"
        "assert HeatExchangerNetwork is hensmith.HeatExchangerNetwork\n"
        "assert HXN2 is hensmith.HeatExchangerNetwork\n"
        "assert ns['HeatExchangerNetwork'] is hensmith.HeatExchangerNetwork\n"
    )

if __name__ == '__main__':
    test_biosteam_first_import_order()
    test_hensmith_first_import_order()
    test_from_imports_and_star_import()
