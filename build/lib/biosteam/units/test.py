# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 17:00:13 2019

@author: Guest Group
"""

from biosteam import Species, Stream, HXprocess
Stream.species = Species('Water', 'Ethanol')

# Simulate counter-current heat exchanger
in1 = Stream('in1', Ethanol=50, T=351.43, phase='l')
in2 = Stream('in2', Water=200, T=373.15, phase='g')
hx = HXprocess('hx', ins=(in1, in2), outs=('out1', 'out2'), Type='ll')
hx.simulate()

# Show results
hx.diagram
hx.show()