# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:13:45 2019

@author: yoelr
"""
import sys
import clr
from System import String
if 'C:\\Users\\yoelr\\AppData\\Local\\DWSIM5' not in sys.path:
    sys.path.append('C:\\Users\\yoelr\\AppData\\Local\\DWSIM5')
clr.AddReference("DWSIM")
clr.AddReference("DWSIM.Thermodynamics.StandaloneLibrary")
from DWSIM.Thermodynamics import StandaloneLibrary as dtl
from DWSIM.Thermodynamics import PropertyPackages
from DWSIM.Thermodynamics import BaseClasses
from DWSIM.Thermodynamics.CalculatorInterface import Calculator
dtlc = Calculator()
dtlc.Initialize()
#dtlc.GetPropPackList()
prpp = dtlc.GetPropPackInstance('Peng-Robinson (PR)')
dtlc.SetupPropertyPackage(prpp, ('Ethanol', 'Water'), (0.5, 0.5))
# myObjectHandle = myAssembly.myClass()
# myObjectHandle.bar()