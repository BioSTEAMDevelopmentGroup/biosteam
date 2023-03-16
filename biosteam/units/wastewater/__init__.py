#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


from ._high_rate_utils import *
from ._high_rate_internal_circulation_rx import *
from ._high_rate_wwt_pump import *
from ._high_rate_polishing_filter import *
from ._high_rate_membrane_bioreactor import *
from ._high_rate_sludge_handling import *

from .conventional import create_conventional_wastewater_treatment_system
from .high_rate import create_high_rate_wastewater_treatment_system

def create_wastewater_treatment_system(WWT='conventional', **WWT_kwargs):
    '''
    Create a wastewater treatment (WWT) system.
    Two configurations are available: conventional and high-rate.
    
    The conventional configuration is based on the WWT process describe in Humbird et al.,
    while the high-rate process is based on the one described in Li et al.
    
    Both configurations are designed to use biological processes to remove 
    organics in the wastewater, and use reverse osmosis to remove the ions, 
    after which the treated water is recycled as process water to be used in the biorefinery.
    
    The main differences between the configurations is that the conventional process
    uses large anaerobic and aerobic basins for organic removal, 
    which has a very low organic loading (grams COD/L/d).
    While the high-rate process uses high-rate anaerobic technologies including
    internal circulation and anaerobic membrane bioreactor with a much higher
    organic loading (>10X of the conventional ones).
    
    With the much higher organic loading rate (smaller reactor and physical footprint)
    and the use of two anaerobic unit operations 
    (anaerobic process does not need aeration therefore saves electricity),
    the high-rate process is expected to have lower CAPEX/OPEX with 
    lower environmental impacts.
    
    Please see the specific functions below for more details on input parameters.
    
    Parameters
    ----------
    WWT : str
        Either "conventional" for the conventional configuration,
        or "high-rate" for the high-rate process.
        All remaining kwargs will be passed to either 
        `create_conventional_wastewater_treatment_system` or
        `create_high_rate_wastewater_treatment_system` depending on the configuration.

    See Also
    --------
    :func:`biosteam.create_conventional_wastewater_treatment_system`;
    :func:`biosteam.create_high_rate_wastewater_treatment_system`;
    '''
    if WWT == 'conventional':
        return create_conventional_wastewater_treatment_system(**WWT_kwargs)
    elif 'high' in WWT and 'rate' in WWT:
        return create_high_rate_wastewater_treatment_system(**WWT_kwargs)
    else:
        raise ValueError(f'The `WWT` value of "{WWT}" is not valid, '
                         'please choose either "conventional" or "high-rate".')
        

from . import (
    _high_rate_utils,
    _high_rate_internal_circulation_rx,
    _high_rate_wwt_pump,
    _high_rate_polishing_filter,
    _high_rate_membrane_bioreactor,
    _high_rate_sludge_handling,
    conventional,
    high_rate,
    )


__all__ = (
    'create_wastewater_treatment_system',
    'create_conventional_wastewater_treatment_system',
    'create_high_rate_wastewater_treatment_system',
    *_high_rate_utils.__all__,
    *_high_rate_internal_circulation_rx.__all__,
    *_high_rate_wwt_pump.__all__,
    *_high_rate_polishing_filter.__all__,
    *_high_rate_membrane_bioreactor.__all__,
    *_high_rate_sludge_handling.__all__,
    # *conventional.__all__,
    # *high_rate.__all__,
    )