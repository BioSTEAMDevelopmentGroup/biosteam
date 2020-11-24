# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 00:16:24 2020

@author: yrc2
"""
from biosteam import units
import biosteam as bst
from biosteam.process_tools import UnitGroup
from thermosteam import functional as fn
from biosteam import main_flowsheet as F
import thermosteam as tmo
import numpy as np

__all__ = ('test_example_sugarcane_subsystem',
)

def test_example_sugarcane_subsystem():
    from biosteam.examples import ethanol_subsystem as ethanol_sys
    biorefinery = UnitGroup('Biorefinery', ethanol_sys.units)
    assert np.allclose(biorefinery.get_installed_cost(), 13.546790896289952, rtol=1e-1)
    assert np.allclose(biorefinery.get_heating_duty(), 149.85496552542645, rtol=1e-2)
    assert np.allclose(biorefinery.get_cooling_duty(), 98.99491946742079, rtol=1e-2)
    assert np.allclose(biorefinery.get_electricity_consumption(), 0.40116044439505544, rtol=1e-2)
    assert np.allclose(biorefinery.get_electricity_production(), 0., rtol=1e-2)
    