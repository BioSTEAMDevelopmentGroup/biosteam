# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 12:37:43 2019

@author: yoelr
"""
import os
import julia
import subprocess
import sys

def install(julia_filepath):
    julia.install(julia=julia_filepath)
    
# try:
#     from julia import Main
# except:
#     filepath = input("Julia installation not found please input file path:")
#     install(filepath)

dirpath = os.path.dirname(os.path.realpath(__file__))
subprocess.check_call(['julia', os.path.join(dirpath, 'vapor_pressure.jl')])
# Main.include(os.path.join(dirpath, "setup.jl"))

# from julia import vapor_pressure
# sys.modules[__name__] = vapor_pressure   # mutate myself