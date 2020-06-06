# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 16:17:00 2017

@author: Yoel Cortes-Pena
"""
from setuptools import setup
#from Cython.Build import cythonize
#import numpy

setup(
    name='biosteam',
    packages=['biosteam'],
    license='MIT',
    version='2.19.0',
    description='The Biorefinery Simulation and Techno-Economic Analysis Modules',
    long_description=open('README.rst').read(),
    author='Yoel Cortes-Pena',
    install_requires=['IPython>=7.9.0', 'biorefineries>=2.14.0',
                      'thermosteam>=0.19.0', 'graphviz>=0.8.3',
                      'chaospy>=3.0.11', 'pipeml>=0.1'],
    python_requires=">=3.6",
    package_data=
        {'biosteam': ['_report/*',
                      '_digraph/*',
                      'utils/*',
                      'compounds/*',
                      'reaction/*',
                      'tests/*',
                      'evaluation/*', 
                      'evaluation/evaluation_tools/*',
                      'process_tools/*',
                      'units/*',
                      'units/design_tools/*',
                      'units/facilities/*',
                      'units/decorators/*',
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