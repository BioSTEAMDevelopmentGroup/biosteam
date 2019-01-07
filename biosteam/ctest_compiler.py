# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 19:28:14 2018

@author: yoelr
"""

from distutils.core import setup
from Cython.Build import cythonize
setup(ext_modules=cythonize('ddl_interop_test.pyx'))