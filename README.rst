========
BioSTEAM
========

.. image:: http://img.shields.io/pypi/v/biosteam.svg?style=flat
   :target: https://pypi.python.org/pypi/biosteam
   :alt: Version_status
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://biosteam.readthedocs.io/en/latest/
   :alt: Documentation
.. image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/yoelcortes/biosteam/blob/master/LICENSE.txt
   :alt: license


.. contents::

What is BioSTEAM?
-----------------

BioSTEAM is an open source process simulation package in Python for fast and flexible preliminary technoeconomic analysis. It was created and developed by `Yoel Cortes-Pena <http://engineeringforsustainability.com/yoelcortespena>`__ and the `Guest Group <http://engineeringforsustainability.com/>`__ as part of the `Center for Advanced Bioenergy and Bioproducts (CABBI) <https://cabbi.bio/>`__ at the `University of Illinois at Urbana-Champaign (UIUC) <https://illinois.edu/>`__. 


Installation
------------

Get the latest version of BioSTEAM from
https://pypi.python.org/pypi/biosteam/

If you have an installation of Python with pip, simple install it with:

    $ pip install biosteam

To get the git version, run:

    $ git clone git://github.com/yoelcortes/biosteam

Documentation
-------------

BioSTEAM's documentation is available on the web:

    http://biosteam.readthedocs.io/

Getting started
---------------

BioSTEAM objects serve as basic building blocks to design and simulate a biorefinery. These include objects that handle material properties, material flows, unit operations, recycle loops and process specifications. To get started, follow the example below:

**Specify the working species** for all streams:

.. code-block:: python
     
   >>> from biosteam import Species, Stream
   >>> Stream.species = Species('Ethanol', 'Water') 

**Create a stream** passing an ID, flow rates, phase, temperature and pressure:

.. code-block:: python

   >>> s1 = Stream(ID='s1', flow=(1, 2), units='kmol/hr', T=300, P=101325)
   >>> s1.show()
   Stream: s1
    phase: 'l', T: 300.00 K, P: 101325 Pa
    flow (kmol/hr): Ethanol   1
                    Water     2
   
**Create a Unit class** passing an ID, IDs for output streams, and any key word arguments specific to the Unit:

.. code-block:: python

   >>> from biosteam.Units import Flash
   >>> # Specify vapor fraction and pressure conditions
   >>> F1 = Flash(ID='F1', outs_ID=('vapor', 'liquid'), V=0.5, P=101325)

**Set the input streams**, and check connections as a Graphviz diagram:

.. code-block:: python

   >>> F1.ins = s1
   >>> F1.diagram

.. figure:: F1_flash.png

**Run the unit** and see the results:

.. code-block:: python

   >>> F1.run()
   >>> F1.show()
   Flash: F1
   ins...
   [0] s1
       phase: 'l', T: 300.00 K, P: 101325 Pa
       flow (kmol/hr): Ethanol   1
                       Water     2
   outs...
   [0] vapor
       phase: 'g', T: 357.20 K, P: 101325 Pa
       flow (kmol/hr): Ethanol   0.761
                       Water     0.739
   [1] liquid
       phase: 'l', T: 357.20 K, P: 101325 Pa
       flow (kmol/hr): Ethanol   0.239
                       Water     1.26

For a more detailed example check out BioSTEAM's documentation.


Bug reports
-----------

To report bugs, please use the BioSTEAM's Bug Tracker at:

    https://github.com/yoelcortes/biosteam


License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the BioSTEAM license, if it is convenient for you,
please cite BioSTEAM if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


Citation
--------

To cite BioSTEAM in publications use::

    Yoel Cortes-Pena (2018). BioSTEAM: The Open-Source Bioprocess Simulation and Technoeconomic Analysis Modules.
    https://github.com/yoelcortes/biosteam
