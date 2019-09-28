# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 23:03:19 2019

This module defines the composition_balance function, which performs mass energy balance the given oil composition of lipid cane to determine its composition.

@author: yoelr
"""
from biosteam import Stream
from biosteam.biorefineries.lipidcane.species import pretreatment_species
from biosteam.biorefineries.lipidcane.system.pretreatment import lipid_cane
from array_collections import tuple_array

__all__ = ('set_lipid_fraction',)

Stream.species = species = pretreatment_species
getattr_ = getattr
carbs_IDs = ('Glucose', 'Sucrose')
fiber_IDs = ('Lignin', 'Cellulose', 'Hemicellulose')
lipid_IDs = ('Lipid',)
water_IDs = ('Water',)

indices = lipid_cane.indices
carbs_index = indices(carbs_IDs)
fiber_index = indices(fiber_IDs)
lipid_index = indices(lipid_IDs)
water_index = indices(water_IDs)

lc = lipid_cane
mol     = lc.mol
carbs   = Stream('Carbs', mol[carbs_index], carbs_IDs)
fiber   = Stream('Fiber', mol[fiber_index], fiber_IDs)
lipid   = Stream('Lipid', mol[lipid_index], lipid_IDs)
streams = (carbs, fiber, lipid)

# Mass property arrays
carbs_mass = lc.mass[carbs_index]
fiber_mass = lc.mass[fiber_index]
lipid_mass = lc.mass[lipid_index]

# Net weight
carbs_massnet = carbs_mass.sum()
fiber_massnet = fiber_mass.sum()
lipid_massnet = lipid_mass.sum()

# Heats of combustion per kg
Hc_carbs_kg = carbs.Hc/carbs_massnet
Hc_lipid_kg = lipid.Hc/lipid_massnet

# Composition
carbs_massfrac = tuple_array(carbs_mass/carbs_massnet)
fiber_massfrac = tuple_array(fiber_mass/fiber_massnet)

water_index, solids_index, ash_index = lc.indices(['Water', 'Solids', 'Ash'])
water_mass = lc.mass.item(water_index)
solids_mass = lc.mass.item(solids_index)
ash_mass = lc.mass.item(ash_index)

def set_lipid_fraction(lipid_fraction):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
    netmass = lc.massnet
    dryfrac = 0.3
    new_lipid_massfrac = lipid_fraction * dryfrac 
    new_carbs_massfrac = 0.149 - new_lipid_massfrac*Hc_lipid_kg/Hc_carbs_kg
    new_fiber_massfrac = dryfrac - new_carbs_massfrac - new_lipid_massfrac - 0.006 - 0.015
    
    water_mass.value = netmass*0.7
    ash_mass.value = netmass*0.006
    solids_mass.value = netmass*0.015
    lipid_mass[:] = new_lipid_massfrac * netmass
    carbs_mass[:] = carbs_massfrac * new_carbs_massfrac * netmass
    fiber_mass[:] = fiber_massfrac * new_fiber_massfrac * netmass
    if any(lc.mol < 0):
        raise ValueError(f'lipid cane oil composition of {lipid_fraction*100:.0f}% dry weight is infeasible')


### Old attempt to include energy balance. It does not work because fiber is more energy dense than sugar (to my surprise) ###

# Hc_0 = lc.Hc
# def Hc_error(carbs_massnet, dryweight_no_oil):
#     """Return the difference between the original heat of combustion and the heat of combustion at the new carbohydrate composition"""
#     carbs_mass[:] = carbs_massfrac * carbs_massnet
#     fiber_mass[:] = fiber_massfrac * (dryweight_no_oil - carbs_massnet)
#     return Hc_0 - lc.Hc

# def composition_balance(oilfrac):
#     """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)."""
#     global carbs_massnet
#     arg = oilfrac
#     oil_mass[:] = mass_oil = oilfrac * dryweight
#     dryweight_no_oil = dryweight - mass_oil
#     carbs_massnet = newton(Hc_error, carbs_massnet, args=(dryweight_no_oil,))
#     if any(lc.mol < 0):
#         raise ValueError(f'Lipid cane oil composition of {arg*100:.0f}% dry weight is infeasible')
#     return carbs_massnet
     
# composition_balance(0.1)
# check = lambda: oil_mass.sum()/dryweight
    