# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Material and specification factors as compiled by [1]_.

References
----------

.. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
    Ng, M. K. Cost Accounting and Capital Cost Estimation.
    In Product and Process Design Principles; Wiley, 2017; pp 426–485.

"""
__all__ = ('vessel_material_factors',
           'pressure_vessel_material_factors',
           'material_densities_lb_per_ft3',
           'material_densities_lb_per_in3',
           'pump_material_factors',
           'pump_gear_factors',
           'pump_centrifugal_factors',
           'distillation_tray_type_factor',
           'tray_material_factor_functions',
           'distillation_column_material_factors',
           'shell_and_tube_material_factor_coefficients',
)

# %% Vessels

#: Material factors for pressure vessels
pressure_vessel_material_factors = {
    'Carbon steel': 1.0,
    'Low-alloy steel': 1.2,
    'Stainless steel 304': 1.7,
    'Stainless steel 316': 2.1,
    'Carpenter 20CB-3': 3.2,
    'Nickel-200': 5.4,
    'Monel-400': 3.6,
    'Inconel-600': 3.9,
    'Incoloy-825': 3.7,
    'Titanium': 7.7
}
#: Material factors for ordinary vessels
vessel_material_factors = {
    'Carbon steel': 1.0,
    'Copper': 1.2,
    'Stainless steel': 2.0,
    'Nickel': 2.5,
    'Titanium clad': 3.0,
    'Titanium': 6.0
}
#: Material densities in lb/ft^3
material_densities_lb_per_ft3 = {
    'Carbon steel': 490,
    'Low-alloy steel': None,
    'Stainless steel 304': 499.4,
    'Stainless steel 316': 499.4,
    'Carpenter 20CB-3': None,
    'Nickel-200': None,
    'Monel-400': None,
    'Inconel-600': None,
    'Incoloy-825': None,
    'Titanium': None
}
#: Material densities in lb/in^3
material_densities_lb_per_in3 = {
    'Carbon steel': 0.284 ,
    'Low-alloy steel': None,
    'Stainless steel 304': 0.289,
    'Stainless steel 316': 0.289,
    'Carpenter 20CB-3': None,
    'Nickel-200': None,
    'Monel-400': None,
    'Inconel-600': None,
    'Incoloy-825': None,
    'Titanium': None
}
# %% Pumps

#: Material factors for pumps
pump_material_factors = {
    'Cast iron':       1,
    'Ductile iron':    1.15,
    'Cast steel':      1.35,
    'Bronze':          1.9,
    'Stainless steel': 2,
    'Hastelloy C':     2.95,
    'Monel':           3.3,
    'Nickel':          3.5,
    'Titanium':        9.7
}
#: Gear-type cost factors for pumps
pump_gear_factors = {
    'OpenDripProof':           1,
    'EnclosedFanCooled':       1.4,
    'ExplosionProofEnclosure': 1.8
}
#: Centrifugal-type cost factors for pumps.
#: Keys are case-split orientation and shaft rpm.
pump_centrifugal_factors = {
    'VSC3600':   1,
    'VSC1800':   1.5,
    'HSC3600':   1.7,
    'HSC1800':   2,
    '2HSC3600':  2.7,
    '2+HSC3600': 8.9
}

# %% Distillation

#: Tray-type cost factors for distillation columns.
distillation_tray_type_factor = {
    'Sieve': 1,
    'Valve': 1.18,
    'Bubble cap': 1.87
}
# Tray material factors (inner diameter, Di, in ft)
def compute_carbon_steel_material_factor(Di):
    return 1

def compute_stainless_steel_304_material_factor(Di):
    return 1.189 + 0.058*Di

def compute_stainless_steel_316_material_factor(Di):
    return 1.401 + 0.073*Di

def compute_carpenter_20CB3_material_factor(Di):
    return 1.525 + 0.079*Di

def compute_monel_material_factor(Di):
    return 2.306 + 0.112*Di

#: Material cost factors for distillation column trays.
tray_material_factor_functions = {
    'Carbon steel': compute_carbon_steel_material_factor,
    'Stainless steel 304': compute_stainless_steel_304_material_factor,
    'Stainless steel 316': compute_stainless_steel_316_material_factor,
    'Carpenter 20CB-3': compute_carpenter_20CB3_material_factor,
    'Monel': compute_monel_material_factor
}
#: Material cost factors for distillation columns.
distillation_column_material_factors = {
    'Carbon steel': 1.0,
    'Low-alloy steel': 1.2,
    'Stainless steel 304': 1.7,
    'Stainless steel 316': 2.1,
    'Carpenter 20CB-3': 3.2,
    'Nickel-200': 5.4,
    'Monel-400': 3.6,
    'Inconel-600': 3.9,
    'Incoloy-825': 3.7,
    'Titanium': 7.7
}

# %% Shell & tube heat exchangers

#: dict[str, tuple[float, float]] Material factor coefficients of shell and 
#: tube heat exchangers. Keys are materials of construction of
#: shell/tube and coefficients are `a` and `b` in Eq. (16.44).
shell_and_tube_material_factor_coefficients =  {
    'Carbon steel/carbon steel':        (0,    0),
    'Carbon steel/brass':	            (1.08, 0.05),
    'Carbon steel/stainles steel':	    (1.75, 0.13),
    'Carbon steel/Monel':	            (2.1,  0.13),
    'Carbon steel/titanium':	        (5.2,  0.16),
    'Carbon steel/Cr-Mo steel':         (1.55, 0.05),
    'Cr-Mo steel/Cr-Mo steel':	        (1.7,  0.07),
    'Stainless steel/stainless steel':  (2.7,  0.07),
    'Monel/Monel':	                    (3.3,  0.08),
    'Titanium/titanium':	            (9.6,  0.06)
}

def compute_shell_and_tube_material_factor(A, a, b):
    r"""
    Return the material factor for shell and tubes given the area [A; ft^3]
    and material-dependent coefficients `a` and `b`. 
    
    Notes
    -----
    Material factors are computed using Eq. 16.44 in [1]_:
    
    :math:`F_M = a + \left(\frac{A}{100.0}\right)^b`
    
    """
    return a + (A/100.0)**b
