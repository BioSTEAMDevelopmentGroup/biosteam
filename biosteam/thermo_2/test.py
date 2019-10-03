# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 12:37:43 2019

@author: yoelr
"""
from biosteam.thermo_2.heat_capacity import HeatCapacityGas, HeatCapacityLiquid, HeatCapacitySolid
from biosteam.thermo_2.vapor_pressure import VaporPressure
from biosteam.thermo_2.volume import VolumeLiquid
from biosteam import Chemical
water = Chemical('Water')
Cpg = HeatCapacityGas(water.CAS, water.MW, water.similarity_variable, False)
Cpl = HeatCapacityLiquid(water.CAS, water.Tc, water.omega,
                         water.MW, water.similarity_variable,
                         water.Cpgm)
Cps = HeatCapacitySolid(water.CAS, water.similarity_variable, water.MW)
VP = VaporPressure(water.CAS)
Vl = VolumeLiquid(water.CAS, Tc=water.Tc, Pc=water.Pc, Tb=water.Tb, Vc=water.Vc, Zc=water.Zc,
                  dipole=water.dipole, Psat=VP)


# import os
# import julia
# import subprocess
# import sys

# def install(julia_filepath):
#     julia.install(julia=julia_filepath)
    
# # try:
# #     from julia import Main
# # except:
# #     filepath = input("Julia installation not found please input file path:")
# #     install(filepath)

# dirpath = os.path.dirname(os.path.realpath(__file__))
# subprocess.check_call(['julia', os.path.join(dirpath, 'vapor_pressure.jl')])
# # Main.include(os.path.join(dirpath, "setup.jl"))

# # from julia import vapor_pressure
# # sys.modules[__name__] = vapor_pressure   # mutate myself