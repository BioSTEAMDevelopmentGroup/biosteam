#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains all wastewater treatment systems and associated unit operations.
The main interface to create a wastewater treatment system is the
:func:`~biosteam.wastewater.create_wastewater_treatment_system` function. Please checkout each submodule for
additional details on the available wastewater treatment systems.

.. toctree::
   :maxdepth: 1
   
   high_rate
   conventional

.. autofunction:: biosteam.wastewater.create_wastewater_treatment_system
   
"""
from .high_rate import create_high_rate_wastewater_treatment_system
from .conventional import create_conventional_wastewater_treatment_system

def create_wastewater_treatment_system(*args, kind=None, **kwargs):
    """
    Create a wastewater treatment system. Two configurations are available: 
    conventional and high-rate. The conventional configuration is based on the 
    process described in Humbird et al., while the high-rate process is based 
    on the one described in Li et al. Both configurations are designed to use
    biological processes to remove organics in the wastewater, and use reverse 
    osmosis to remove the ions, after which the treated water is recycled as 
    process water to be used in the biorefinery.
    
    The main differences between the configurations is that the conventional process
    uses large anaerobic and aerobic basins for organic removal, 
    which has a very low organic loading (grams COD/L/d).
    While the high-rate process uses high-rate anaerobic technologies including
    internal circulation and anaerobic membrane bioreactor with a much higher
    organic loading (>10X of the conventional ones). With the much higher organic 
    loading rate (smaller reactor and physical footprint)
    and the use of two anaerobic unit operations 
    (anaerobic process does not need aeration therefore saves electricity),
    the high-rate process is expected to have lower CAPEX/OPEX with 
    lower environmental impacts
    (see the example for MESP comparison for the corn stover biorefinery).
    
    Parameters
    ----------
    *args: 
        These parameters will be passed to either 
        :func:`~biosteam.wastewater.conventional.create_conventional_wastewater_treatment_system` or
        :func:`~biosteam.wastewater.high_rate.system.create_high_rate_wastewater_treatment_system`
        depending on the configuration.
    kind : str, optional
        Either 'conventional' for the conventional configuration,
        or "high-rate" for the high-rate process. Defaults to 'conventional'.
    
    Other Parameters
    ----------------
    **kwargs : 
        All remaining parameters will be passed to either 
        :func:`~biosteam.wastewater.conventional.create_conventional_wastewater_treatment_system` or
        :func:`~biosteam.wastewater.high_rate.system.create_high_rate_wastewater_treatment_system`
        depending on the configuration.
        
    Examples
    --------
    Calculate the minimum ethanol selling price of the corn stover biorefinery 
    with either the conventional or high-rate wastewater treatment configurations:
    
    >>> import biosteam as bst
    >>> from biorefineries import cornstover as cs
    >>> bst.settings.set_thermo(cs.create_chemicals())
    >>> def get_MESP(**WWT_kwargs):
    ...     bst.main_flowsheet.set_flowsheet(WWT_kwargs['kind'])
    ...     # WWT_kwargs will be passed to the wastewater treatment system creator
    ...     sys = cs.create_system(WWT_kwargs=WWT_kwargs) 
    ...     sys.simulate()
    ...     tea = cs.create_tea(sys)
    ...     ethanol = sys.get_outlet('ethanol')
    ...     MESP = tea.solve_price(ethanol) * cs.ethanol_density_kggal
    ...     print(f"{WWT_kwargs['kind']} MESP: ${round(MESP, 1)}/gal")
    
    >>> # With the conventional WWT process
    >>> get_MESP(kind='conventional')
    conventional MESP: $2.0/gal
    
    >>> # With the high-rate WWT process
    >>> get_MESP(process_ID=6, kind='high-rate')
    high-rate MESP: $1.4/gal
    
    References
    ----------
    [1] Li et al., Design of a High-Rate Wastewater Treatment Process 
    for Energy and Water Recovery at Biorefineries. 
    ACS Sustainable Chem. Eng. 2023, 11 (9), 3861–3872. 
    https://doi.org/10.1021/acssuschemeng.2c07139.

    [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

    [3] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234
    
    """
    kind = 'conventional' if kind is None else kind.translate({ord(i): None for i in '_- '}).lower()
    if kind == 'conventional':
        return create_conventional_wastewater_treatment_system(*args, **kwargs)
    elif kind == 'highrate':
        return create_high_rate_wastewater_treatment_system(*args, **kwargs)
    else:
        raise ValueError(f"invalid `kind` '{kind}'; `kind` must be either "
                         "'conventional' or 'high-rate'")

from . import (
    conventional,
    high_rate,
)

__all__ = (
    'create_wastewater_treatment_system',
    'create_conventional_wastewater_treatment_system',
    'create_high_rate_wastewater_treatment_system',
)
