=========================================================================
BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
=========================================================================

.. image:: http://img.shields.io/pypi/v/biosteam.svg?style=flat
   :target: https://pypi.python.org/pypi/biosteam
   :alt: Version_status
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://biosteam.readthedocs.io/en/latest/
   :alt: Documentation
.. image:: http://img.shields.io/badge/license-UIUC-blue.svg?style=flat
   :target: https://github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
   :alt: license
.. image:: https://img.shields.io/pypi/pyversions/biosteam.svg
   :target: https://pypi.python.org/pypi/biosteam
   :alt: Supported_versions
.. image:: https://zenodo.org/badge/164639830.svg
   :target: https://zenodo.org/badge/latestdoi/164639830
.. image:: https://coveralls.io/repos/github/BioSTEAMDevelopmentGroup/biosteam/badge.svg?branch=master
   :target: https://coveralls.io/github/BioSTEAMDevelopmentGroup/biosteam?branch=master
.. image:: https://badges.gitter.im/BioSTEAM-users/BioSTEAM.svg
   :alt: Join the chat at https://gitter.im/BioSTEAM-users/community
   :target: https://gitter.im/BioSTEAM-users/community

.. contents::

What is BioSTEAM?
-----------------

BioSTEAM is a fast and flexible package for the design, simulation, 
techno-economic analysis, and life cycle assessment of biorefineries under uncertainty [1]_. 
BioSTEAM is built to streamline and automate early-stage technology evaluations 
and to enable rigorous sensitivity and uncertainty analyses. Complete 
biorefinery configurations are available at the `Bioindustrial-Park 
<https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`_ GitHub repository, 
BioSTEAM's premier repository for biorefinery models and results. The long-term 
growth and maintenance of BioSTEAM is supported through both community-led 
development and the research institutions invested in BioSTEAM, including the 
`Center for Advanced Bioenergy and Bioproducts Innovation (CABBI) <https://cabbi.bio/>`_. 
Through its open-source and community-lead platform, BioSTEAM aims to foster 
communication and transparency within the biorefinery research community for an 
integrated effort to expedite the evaluation of candidate biofuels and 
bioproducts.

All data on chemicals, fitted coefficients, and functions to estimate chemical 
properties are imported from `chemicals <https://github.com/CalebBell/chemicals>`_,
a community-driven open-source library developed by Caleb Bell. The 
`thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_ library, 
BioSTEAM's premire thermodynamic engine, creates thermodynamic property packages 
by implementing both the mixing rules and the algorithms for thermodynamic equilibrium 
and by relying on the `chemicals` library for pure component models.

Installation
------------

Get the latest version of BioSTEAM from `PyPI <https://pypi.python.org/pypi/biosteam/>`__. If you have an installation of Python with pip, simple install it with:

    $ pip install biosteam

To get the git version, run:

    $ git clone git://github.com/BioSTEAMDevelopmentGroup/biosteam

For help on common installation issues, please visit the `documentation <https://biosteam.readthedocs.io/en/latest/#installation>`__.

Documentation
-------------

BioSTEAM's documentation is available on the web:

    http://biosteam.readthedocs.io/

Bug reports
-----------

To report bugs, please use the BioSTEAM's Bug Tracker at:

    https://github.com/BioSTEAMDevelopmentGroup/biosteam

Contributing
------------
For guidelines on how to contribute visit:

    https://biosteam.readthedocs.io/en/latest/contributing/index.html

We also welcome anyone with an interest in discussing contributions to BioSTEAM to join our
monthly meeting at 8-9:30am CDT, every 3rd Friday. Please email biosteamdevelopmentgroup@gmail.com 
for a link to join.

License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the BioSTEAM license, if it is convenient for you,
please cite BioSTEAM if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


About the authors
-----------------

BioSTEAM was created and developed by `Yoel Cortes-Pena <http://engineeringforsustainability.com/yoelcortespena>`__ as part of the `Guest Group <http://engineeringforsustainability.com/>`__ and the `Center for Advanced Bioenergy and Bioproducts Innovation (CABBI) <https://cabbi.bio/>`__ at the `University of Illinois at Urbana-Champaign (UIUC) <https://illinois.edu/>`__. 

References
----------
.. [1] `Cortes-Peña et al. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020. <https://doi.org/10.1021/acssuschemeng.9b07040>`__.


