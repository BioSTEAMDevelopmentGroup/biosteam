# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:17:00 2017

@author: Yoel Cortes-Pena
"""
from setuptools import setup

setup(
    name='biosteam',
    packages=['biosteam'],
    license='MIT',
    version='0.01',
    description='The Open-Source Bioprocess Simulation and Techno-Economic Analysis Modules',
    long_description=open('README.rst').read(),
    author='Yoel Cortes-Pena',
    install_requires=['pint', 'scipy', 'IPython', 'thermo', 'bookkeep', 'colorpalette',
                      'pandas', 'cashflows', 'numpy', 'graphviz', 'matplotlib', 'openpyxl'],
    package_data={'biosteam': ['my_units_defs.txt']},
    platforms=["Windows"],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/yoelcortes/biosteam',
    download_url='https://github.com/yoelcortes/biosteam.git',
)

