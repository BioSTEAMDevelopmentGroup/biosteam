# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from setuptools import setup
#from Cython.Build import cythonize
#import numpy

setup(
    name='biosteam',
    packages=['biosteam'],
    license='MIT',
    version='2.20.16',
    description='The Biorefinery Simulation and Techno-Economic Analysis Modules',
    long_description=open('README.rst').read(),
    author='Yoel Cortes-Pena',
    install_requires=['IPython>=7.9.0', 'biorefineries>=2.15.11',
                      'thermosteam>=0.20.11', 'graphviz>=0.8.3',
                      'chaospy>=3.0.11', 'pipeml>=0.1'],
    python_requires=">=3.6",
    package_data=
        {'biosteam': ['report/*',
                      'digraph/*',
                      'utils/*',
                      'compounds/*',
                      'reaction/*',
                      'tests/*',
                      'evaluation/*', 
                      'evaluation/evaluation_tools/*',
                      'process_tools/*',
                      'plots/*',
                      'units/*',
                      'units/design_tools/*',
                      'units/facilities/*',
                      'units/decorators/*',
                      'examples/*',
                      ]},
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/BioSTEAMDevelopmentGroup/biosteam',
    download_url='https://github.com/BioSTEAMDevelopmentGroup/biosteam.git',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'License :: OSI Approved :: University of Illinois/NCSA Open Source License',
                 'Programming Language :: Python :: 3.6',
				 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics'],
    keywords='chemical process simmulation bioprocess engineering mass energy balance material properties phase equilibrium CABBI biorefinery biofuel bioproducts',
)