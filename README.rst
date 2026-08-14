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

**Read in:** `Español <README.es.rst>`_

.. contents::

What is BioSTEAM?
-----------------

BioSTEAM is a fast and flexible package for the design, simulation, 
techno-economic analysis, and life cycle assessment of chemical processes under uncertainty [1]_. 
BioSTEAM is built to streamline and automate early-stage technology evaluations 
and to enable rigorous sensitivity and uncertainty analyses. BioSTEAM has been 
used to model a range of fermentation-based bioproduct pathways
going from conventional/cellulosic crops, municipal solid waste, and flue gas all
the way to diols, organic acids, oleochemicals, and biofuels. Complete 
biorefinery configurations are available at the `Bioindustrial-Park 
<https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`_ GitHub repository, 
BioSTEAM's premier repository for biorefinery models and results. Newer applications
for BioSTEAM include thermochemical upcycling of waste plastics through
pyrolysis and solvent-based dissolution and precipitation. 

The long-term growth and maintenance of BioSTEAM is supported through both community-led 
development and the research institutions invested in BioSTEAM. 
Through its open-source and community-lead platform, BioSTEAM aims to foster 
communication and transparency within the process systems research community for an 
integrated effort to expedite the development of sustainable processes.

Data on chemicals and algorithms to estimate thermodynamic properties are 
imported from `chemicals <https://github.com/CalebBell/chemicals>`__
and `thermo <https://github.com/CalebBell/chemicals>`__,
community-driven open-source libraries developed by Caleb Bell. BioSTEAM's 
premire thermodynamic engine, `ThermoSTEAM <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`__, 
builds upon these libraries to facilitate the creation of thermodynamic property packages.

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
For guidelines on how to contribute, visit:

    https://biosteam.readthedocs.io/en/latest/contributing/index.html


License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the BioSTEAM license, if it is convenient for you,
please cite BioSTEAM if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


About the authors
-----------------

BioSTEAM was originally created and developed by `Yoel Cortés-Peña <https://blog.nus.edu.sg/corteslab/>`__ 
as part of the `Guest Group <http://engineeringforsustainability.com/yoelcortespena>`__ 
and the `Center for Advanced Bioenergy and Bioproducts Innovation (CABBI) <https://cabbi.bio/>`__ 
at the `University of Illinois at Urbana-Champaign (UIUC) <https://illinois.edu/>`__. 
Yoel is now an Assistant Professor at the National University of Singapore,
where he and his `lab group <https://blog.nus.edu.sg/corteslab/>`__ 
continue to develop BioSTEAM's core simulation capabilities. 

BioSTEAM features contributons from numerous other lab groups and 
individuals. Here is the list of other core contributers who have committed significant
resources, expertise, and other impactful efforts:

* `Jeremy Guest <http://engineeringforsustainability.com/people/>`__ led efforts to found, develop, coordinate, and expand the software.

* `Sarang Sunil Bhagwat <https://github.com/sarangbhagwat>`__ led the development of the heat exchanger network and contributed enhancements to unit operations.

* `Yalin Li <https://yalinli.group/>`__ contributed wastewater treatment models and general enhancements.

* `Rui Shi <https://www.linkedin.com/in/chuyingshi/>`__ contributed kinetic fermentation model.

BioSTEAM has also received numerous contributions from the community. You can view
direct contributions and developers through this 
`GitHub link <https://github.com/BioSTEAMDevelopmentGroup/biosteam/graphs/contributors>`__.

References
----------
.. [1] `Cortés-Peña et al. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020. <https://doi.org/10.1021/acssuschemeng.9b07040>`__.


