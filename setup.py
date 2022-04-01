# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
from setuptools import setup

setup(
    name='biosteam',
    packages=['biosteam'],
    license='MIT',
    version='2.32.4',
    description='The Biorefinery Simulation and Techno-Economic Analysis Modules',
    long_description=open('README.rst', encoding='utf-8').read(),
    author='Yoel Cortes-Pena',
    install_requires=['IPython>=7.9.0',
                      'thermosteam>=0.28.6', 
                      'graphviz>=0.17',
                      'chaospy>=3.3.9'],
    extras_require={ 
        'dev': [
            'biorefineries>=2.23.10',
            'sympy',
            'sphinx', 
            'sphinx_rtd_theme', 
            'pyyaml',
            'pytest-cov',
            'coveralls',
            'SALib',
        ]
    }, 
    package_data=
        {'biosteam': ['report/*',
                      'digraph/*',
                      'utils/*',
                      'evaluation/*', 
                      'evaluation/evaluation_tools/*',
                      'process_tools/*',
                      'plots/*',
                      'units/*',
                      'units/facilities/hxn/*',
                      'units/design_tools/*',
                      'units/facilities/*',
                      'units/decorators/*',
                      ]},
    python_requires='>=3.8',
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/BioSTEAMDevelopmentGroup/biosteam',
    download_url='https://github.com/BioSTEAMDevelopmentGroup/biosteam.git',
    classifiers=['License :: OSI Approved :: University of Illinois/NCSA Open Source License',
                 'Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Intended Audience :: Developers',
                 'Intended Audience :: Education',
                 'Intended Audience :: Manufacturing',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: POSIX :: BSD',
                 'Operating System :: POSIX :: Linux',
                 'Operating System :: Unix',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: Implementation :: CPython',
                 'Topic :: Education'],
    keywords=['chemical process simmulation', 'bioprocess engineering', 'mass and energy balance', 'material properties', 'phase equilibrium', 'CABBI', 'biorefinery', 'biofuel', 'bioproducts'],
)