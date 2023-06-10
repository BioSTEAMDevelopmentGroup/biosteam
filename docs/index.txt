The Biorefinery Simulation and TEA Modules
==========================================

.. toctree::
   :maxdepth: 2
   :hidden:
   
   tutorial/index
   API/index
   contributing/index

.. grid:: 1 1 2 2

    .. grid-item::
    
        .. image:: images/demo_gif_dark.gif
           :class: only-dark
           :align: center

        .. image:: images/demo_gif_light.gif
           :class: only-light
           :align: center


    .. grid-item::

        BioSTEAM is an open-source platform that streamlines the design, simulation, techno-economic 
        analysis (TEA) and life-cycle assessment (LCA) of biorefineries across thousands 
        of scenarios. BioSTEAM is also leveraged by `QSDsan <https://qsdsan.readthedocs.io/en/latest/>`_,
        a library for the quantitative sustainable design of sanitation and resource recovery
        systems. The long-term growth and maintenance of BioSTEAM is supported
        through both community-led development and the research institutions invested in BioSTEAM, 
        including the `Center for Advanced Bioenergy and Bioproducts Innovation <https://cabbi.bio/>`_.


.. grid:: 1 2 3 4
   
   
    .. grid-item-card:: Getting Started
       :text-align: center
       :link: https://biosteam.readthedocs.io/en/latest/tutorial/Getting_started.html
       :link-type: url
       :padding: 1

       .. figure:: images/aerobic_chamber.png
          :height: 100
          :class: dark-light
          :align: center
          
       Tutorial on BioSTEAM's features.


    .. grid-item-card:: Bioindustrial-Park
       :text-align: center
       :link: https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park
       :link-type: url
       :padding: 1
       
       .. figure:: images/boiler.png
          :height: 100
          :class: dark-light
          :align: center
       
       Repository for biorefinery models.

        
    .. grid-item-card:: API Reference
       :text-align: center
       :link: https://biosteam.readthedocs.io/en/latest/API/index.html
       :link-type: url
       :padding: 1
       
       .. figure:: images/fermenter.png
          :height: 100
          :class: dark-light
          :align: center
       
       Documentation on all functionality. 
       
       
    .. grid-item-card:: QSDsan
       :text-align: center
       :link: https://qsdsan.readthedocs.io/en/latest/
       :link-type: url
       :padding: 1
       
       .. figure:: images/GG_logo.png
          :height: 100
          :class: dark-light
          :align: center
       
       Sanitation and resource recovery. 

Installation
------------

If you have an installation of Python with pip, simple install it with:

.. code-block:: bash

   $ pip install biosteam


To get the git version, run:

.. code-block:: bash
   
   $ git clone --depth 10 git://github.com/BioSTEAMDevelopmentGroup/biosteam


Common Issues
-------------

* **Cannot install/update BioSTEAM:**

  If you are having trouble installing or updating BioSTEAM, it may be due to dependency issues. You can bypass these using:
  
  .. code-block:: bash

     $ pip install --user --ignore-installed biosteam

  You can make sure you install the right version by including the version number:

  .. code-block:: bash

     $ pip install biosteam==<version>

* **Unit and system diagrams are not displaying:**

  BioSTEAM uses `Graphviz <http://www.graphviz.org/>`__ to make flowsheet diagrams. To properly install Graphviz in an anaconda distribution, please run the following line:
  
  .. code-block:: bash

     $ conda install graphviz

  Additionally, please follow the following instructions for `installing graphviz on windows <https://forum.graphviz.org/t/new-simplified-installation-procedure-on-windows/224#format-svg-not-recognized-use-one-of>`__.

Scientific Papers
-----------------

Several studies have leveraged the BioSTEAM platform to compare conversion 
technologies, chart development pathways for various bioproducts, and build new
tools for sustainability assessments. Here is a short list of related publications:

* **Software tools**:

  - `BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020 <https://doi.org/10.1021/acssuschemeng.9b07040>`_

  - `Thermosteam: BioSTEAM's Premier Thermodynamic Engine. Journal of Open Source Software. 2020 <doi.org/10.21105/joss.02814>`_

  - `QSDsan: an integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. 2022 <https://doi.org/10.1039/D2EW00455K>`_

* **Social, economic, and policy studies**:

  - `An agent-based modeling tool supporting bioenergy and bio-product community communication regarding cellulosic bioeconomy development. Renewable and Sustainable Energy Reviews 2022 <https://doi.org/10.1016/j.rser.2022.112745>`_

  - `Implications of Biorefinery Policy Incentives and Location-Specific Economic Parameters for the Financial Viability of Biofuels. Environ. Sci. Technol. 2023 <https://doi.org/10.1021/acs.est.2c07936>`_

* **Bioproduct and biofuel studies**:

  - `Metabolic engineering of Yarrowia lipolytica to produce 1,2-diacyl-3-acetyl triacylglycerol. Metabolic Engineering. 2023 <https://doi.org/10.1016/j.ymben.2023.01.003>`_

  - `Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries. ACS Sustainable Chem. Eng. 2022 <https://doi.org/10.1021/acssuschemeng.2c04204>`_

  - `Adsorptive separation and recovery of triacetic acid lactone from fermentation broth. Biofuels Bioprod Bioref 2022 <https://doi.org/10.1002/bbb.2427>`_

  - `Rewiring yeast metabolism for producing 2,3-butanediol and two downstream applications: Techno-economic analysis and life cycle assessment of methyl ethyl ketone (MEK) and agricultural biostimulant production. Chemical Engineering Journal 2022 <https://doi.org/10.1016/j.cej.2022.138886>`_

  - `Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021 <https://doi.org/10.1021/acssuschemeng.1c05441>`_
  
  - `Renewable linear alpha-olefins by base-catalyzed dehydration of biologically-derived fatty alcohols. Green Chemistry 2021 <https://doi.org/10.1039/D1GC00243K>`_

  - `Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021 <https://doi.org/10.1021/acssuschemeng.0c08055>`_

  - `Techno-Economic Evaluation of Biorefineries Based on Low-Value Feedstocks Using the BioSTEAM Software: A Case Study for Animal Bedding. Processes 2020 <https://doi.org/10.3390/pr8080904>`_

* **Plastics and recycling**:

  - `High-purity polypropylene from disposable face masks via solvent-targeted recovery and precipitation. Green Chemistry. 2023 <https://doi.org/10.1039/D3GC00205E>`_


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`