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
    version='0.2.1',
    description='The Techno-Economic Analysis Modules',
    long_description=open('README.rst').read(),
    #ext_modules=cythonize('biosteam/equilibrium/unifac.pyx'),
    #include_dirs=[numpy.get_include()],
    author='Yoel Cortes-Pena',
    install_requires=['pint==0.9', 'ht==0.1.52', 'fluids==0.1.74',
                      'scipy', 'IPython', 'thermo==0.1.39', 'colorpalette==0.1.1',
                      'array_collections==0.1.7', 'free_properties==0.1.8',
                      'pandas', 'numpy', 'graphviz', 'matplotlib'],
    python_requires=">=3.6",
    package_data=
        {'biosteam': ['equilibrium/*', 'inspect/*', 'price/*', 'report/*',
                      'my_units_defs.txt', 'utils/*', 'units/*', 'evaluation/*',
                      'units/designtools/*'], },
    platforms=["Windows", "Mac"],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/yoelcortes/biosteam',
    download_url='https://github.com/yoelcortes/biosteam.git',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics'],
    keywords='chemical process simmulation bioprocess engineering mass energy balance material properties phase equilibrium CABBI',
)