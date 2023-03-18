# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Unit construction and functions for creating a high-rate wastewater treatment system
as described in Li et al., with some unit operations based on Humbird et al., and Davis et al.

.. contents:: :local:

Systems
-------
.. autofunction:: biosteam.wastewater.high_rate.system.create_high_rate_wastewater_treatment_system
    
Unit operations
---------------
.. autoclass:: biosteam.wastewater.high_rate.internal_circulation_rx.InternalCirculationRx 
.. autoclass:: biosteam.wastewater.high_rate.membrane_bioreactor.AnMBR
.. autoclass:: biosteam.wastewater.high_rate.polishing_filter.PolishingFilter
.. autoclass:: biosteam.wastewater.high_rate.sludge_handling.BeltThickener
.. autoclass:: biosteam.wastewater.high_rate.sludge_handling.SludgeCentrifuge
.. autoclass:: biosteam.wastewater.high_rate.BiogasUpgrading
.. autoclass:: biosteam.wastewater.high_rate.CHP
.. autoclass:: biosteam.wastewater.high_rate.ReverseOsmosis
.. autoclass:: biosteam.wastewater.high_rate.Skipped

References
----------
.. [1] Li et al., Design of a High-Rate Wastewater Treatment Process 
    for Energy and Water Recovery at Biorefineries. 
    ACS Sustainable Chem. Eng. 2023, 11 (9), 3861–3872. 
    https://doi.org/10.1021/acssuschemeng.2c07139.

.. [2] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

.. [3] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234

.. [4] Kontos, G. A. Advanced Anaerobic Treatment for Energy Recovery and Improved Process Economics
    in the Management of Biorefinery Wastewaters.
    M.S. Thesis, University of Illinois Urbana-Champaign, Urbana, IL, 2021.

.. [5] Schueller, D. MUNICIPAL RESIDENTIAL WASTEWATER RATES.
    Metropolitan Council Environmental Services, 2020.

.. [6] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
    Valorization of Dilute Organic Carbon Waste Streams.
    Energy Environ. Sci. 2016, 9 (3), 1102–1112.
    https://doi.org/10.1039/C5EE03715H.

.. [7] Irizar et al., Model-Based Design of a Software Sensor for Real-Time
    Diagnosis of the Stability Conditions in High-Rate Anaerobic Reactors –
    Full-Scale Application to Internal Circulation Technology.
    Water Research 2018, 143, 479–491.
    `<https://doi.org/10.1016/j.watres.2018.06.055>`_.
    
.. [8] Tchobanoglous et al., Wastewater Engineering: Treatment and Resource Recovery,
    5th ed.; McGraw-Hill Education: New York, 2013.

.. [9] `Industrial filtering equipment gravity thickener rotary thickening belt filter press \
    <https://www.alibaba.com/product-detail/Industrial-filtering-equipment-gravity-thickener-rotary_60757627922.html?spm=a2700.galleryofferlist.normal_offer.d_title.78556be9t8szku>`_
    Data obtained on 7/21/2021.

"""

from .utils import *
from .internal_circulation_rx import *
from .wwt_pump import *
from .polishing_filter import *
from .membrane_bioreactor import *
from .sludge_handling import *
from .system import *

from . import (
    utils,
    internal_circulation_rx,
    wwt_pump,
    polishing_filter,
    membrane_bioreactor,
    sludge_handling,
    system,
)

__all__ = (
    *utils.__all__,
    *internal_circulation_rx.__all__,
    *wwt_pump.__all__,
    *polishing_filter.__all__,
    *membrane_bioreactor.__all__,
    *sludge_handling.__all__,
    *system.__all__,
)
