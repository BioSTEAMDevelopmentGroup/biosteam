# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 00:22:07 2019

@author: yoelr
"""

from biosteam.units import Mixer, Splitter
M1 = Mixer('M1')
S1 = Splitter('S1', split=0.5) # Split to 0th output stream
print(S1)