# -*- coding: utf-8 -*-
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
           'distillation_tray_type_factor')

# %% Vessels

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
    'Titanium': 7.7}

vessel_material_factors = {
    'Carbon steel': 1.0,
    'Copper': 1.2,
    'Stainless steel': 2.0,
    'Nickel': 2.5,
    'Titanium clad': 3.0,
    'Titanium': 6.0}

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
    'Titanium': None}

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
    'Titanium': None}

# %% Pumps

# Material factors
pump_material_factors = {
    'Cast iron':       1,
    'Ductile iron':    1.15,
    'Cast steel':      1.35,
    'Bronze':          1.9,
    'Stainless steel': 2,
    'Hastelloy C':     2.95,
    'Monel':           3.3,
    'Nickel':          3.5,
    'Titanium':        9.7}

# Gear factors
pump_gear_factors = {
    'OpenDripProof':           1,
    'EnclosedFanCooled':       1.4,
    'ExplosionProofEnclosure': 1.8}

pump_centrifugal_factors = {
    'VSC3600':   1,
    'VSC1800':   1.5,
    'HSC3600':   1.7,
    'HSC1800':   2,
    '2HSC3600':  2.7,
    '2+HSC3600': 8.9}


# %% Distillation

distillation_tray_type_factor = {
    'Sieve': 1,
    'Valve': 1.18,
    'Bubble cap': 1.87}

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

tray_material_factor_functions = {
    'Carbon steel': compute_carbon_steel_material_factor,
    'Stainless steel 304': compute_stainless_steel_304_material_factor,
    'Stainless steel 316': compute_stainless_steel_316_material_factor,
    'Carpenter 20CB-3': compute_carpenter_20CB3_material_factor,
    'Monel': compute_monel_material_factor}

# Column Material
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
    'Titanium': 7.7}

# %% Shell & tube heat exchangers

# Materials of Construction Shell/Tube (a and b in Eq. (16.44))
shell_and_tube_material_factor_coefficients =  {
    'Carbon steel/carbon steel':       (0, 0),
    'Carbon steel/brass':	            (1.08, 0.05),
    'Carbon steel/stainles steel':	  (1.75, 0.13),
    'Carbon steel/Monel':	            (2.1, 0.13),
    'Carbon steel/titanium':	       (5.2, 0.16),
    'Carbon steel/Cr-Mo steel':        (1.55, 0.05),
    'Cr-Mo steel/Cr-Mo steel':	       (1.7, 0.07),
    'Stainless steel/stainless steel': (2.7, 0.07),
    'Monel/Monel':	                 (3.3, 0.08),
    'Titanium/titanium':	            (9.6, 0.06)}

def compute_shell_and_tube_material_factor(A, shell_and_tube_material):
    r"""
    Return the material factor for shell and tubes given the area [A; ft^3]
    and the shell and tube material. 
    
    Notes
    -----
    Material factors are computed using Eq. 16.44 in [1]_, which relies on 
    material-dependent coefficients `a` and `b`:
    
    :math:`F_M = a + \left(\frac{A}{100.0}\right)^b`
    
    
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.;
        Ng, M. K. Cost Accounting and Capital Cost Estimation.
        In Product and Process Design Principles; Wiley, 2017; pp 426–485.
    """   
    a, b = shell_and_tube_material_factor_coefficients[shell_and_tube_material]
    return a + (A/100.0)**b
