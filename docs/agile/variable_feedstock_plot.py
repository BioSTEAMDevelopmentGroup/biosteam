# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 03:52:42 2021

@author: yrc2
"""

import matplotlib.pyplot as plt
from biosteam.utils import colors as c

# x-axis; 'Operating time [hr]'
# y-axis; 'Flow rate [kg/s]'
x = [0,    1000, 2000, 3000, 4000]
y = [50,   100,  150,  100,    50]
index = list(range(len(y)))
plt.figure()
plt.bar(index, y, width=1., align='edge', edgecolor=c.CABBI_black.RGBn, color=c.CABBI_blue.RGBn)
plt.xticks(index + [index[-1] + 1], x + [5000])
plt.xlabel('Operating time [hr]')
plt.ylabel('Feedstock flow rate [kg/s]')
plt.show()


# x-axis; 'Operating time [hr]'
# y-axis; 'Flow rate [kg/s]'
x = [0,    500, 1000, 1500, 2000]
y = [20,   40,  60,  40,    20]
index = list(range(len(y)))
plt.figure()
plt.bar(index, y, width=1., align='edge', edgecolor=c.CABBI_black.RGBn, color=c.CABBI_orange.RGBn)
plt.xticks(index + [index[-1] + 1], x + [2500])
plt.xlabel('Operating time [hr]')
plt.ylabel('Feedstock flow rate [kg/s]')
plt.show()