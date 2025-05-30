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

       .. image:: images/phenomenode.png
          :height: 100
          :class: only-light
          :align: center
          
       .. image:: images/phenomenode_dark_2.png
          :height: 100
          :class: only-dark
          :align: center
          
       Tutorials on BioSTEAM


    .. grid-item-card:: Bioindustrial-Park
       :text-align: center
       :link: https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park
       :link-type: url
       :padding: 1
       
       .. image:: images/process.png
          :height: 100
          :class: only-light
          :align: center
          
       .. image:: images/process_dark.png
          :height: 100
          :class: only-dark
          :align: center
       
       Biorefinery models

        
    .. grid-item-card:: API Reference
       :text-align: center
       :link: https://biosteam.readthedocs.io/en/latest/API/index.html
       :link-type: url
       :padding: 1
       
       .. image:: images/fermenter.png
          :height: 100
          :class: dark-light
          :align: center
       
       Detailed documentation
       
    .. grid-item-card:: QSDsan
       :text-align: center
       :link: https://qsdsan.readthedocs.io/en/latest/
       :link-type: url
       :padding: 1
    
       .. image:: images/GG_logo.png
          :height: 100
          :class: dark-light
          :align: center
    
       Sanitation & resources


Installation
------------

#. If you have an installation of Python with pip, simple install it with:

   .. code-block:: bash

      $ pip install biosteam


   To get the git version, run:

   .. code-block:: bash
   
      $ git clone --depth 10 git://github.com/BioSTEAMDevelopmentGroup/biosteam

#. BioSTEAM uses `Graphviz <http://www.graphviz.org/>`__ to make flowsheet diagrams. 
   You will need to install Graphviz separately as follows:

   * Windows: Download the EXE installer and follow the instructions listed `in this link <https://graphviz.org/download/>`__

   * Ubuntu: 

     .. code-block:: bash
    
       $ sudo apt-get install graphviz
   
   * MacOS: 

     .. code-block:: bash
    
        $ brew install graphviz
   
#. To properly install Graphviz in an anaconda distribution, run the following line:

   .. code-block:: bash
    
      $ conda install python-graphviz

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

  Graphviz may not be properly installed or may be missing from your python path. 
  Please follow the graphviz installation procedure outlined above.

Scientific Papers
-----------------

Several studies have leveraged the BioSTEAM platform to compare conversion 
technologies, chart development pathways for various bioproducts, and build new
tools for sustainability assessments. Here is a short list of related publications:

* **Software tools**:

  #. `BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem Eng 2020 <https://doi.org/10.1021/acssuschemeng.9b07040>`_

  #. `Thermosteam: BioSTEAM's Premier Thermodynamic Engine. Journal of Open Source Software 2020 <doi.org/10.21105/joss.02814>`_

  #. `QSDsan: an integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ Sci: Water Res Technol 2022 <https://doi.org/10.1039/D2EW00455K>`_

* **Social, economic, and policy studies**:

  #. `An agent-based modeling tool supporting bioenergy and bio-product community communication regarding cellulosic bioeconomy development. Renewable and Sustainable Energy Reviews 2022 <https://doi.org/10.1016/j.rser.2022.112745>`_

  #. `Implications of Biorefinery Policy Incentives and Location-Specific Economic Parameters for the Financial Viability of Biofuels. Environ Sci Technol 2023 <https://doi.org/10.1021/acs.est.2c07936>`_

* **Bioproduct and biofuel studies**:
  
  #. `Sustainable potassium sorbate production from triacetic acid lactone in food-grade solvents. Green Chem 2025 <https://doi.org/10.1039/D4GC04832F>`_
  
  #. `Integration of plant and microbial oil processing at oilcane biorefineries for more sustainable biofuel production. GCB Bioenergy 2024 <https://doi.org/10.1111/gcbb.13183>`_

  #. `Economic and Environmental Sustainability of Bio-Based HMF Production and Recovery from Lignocellulosic Biomass. Green Chem 2024 <https://doi.org/10.1039/D4GC04270K>`_

  #. `Comparative techno-economic and life cycle assessment of electrocatalytic processes for lignin valorization. Green Chem 2024 <https://doi.org/10.1039/D4GC01963F>`_

  #. `The Economic and Environmental Case for Cattle Manure and Prairie Grass-Derived Sustainable Aviation Fuel. Energy Fuels 2024 <https://doi.org/10.1021/acs.energyfuels.4c02929>`_

  #. `A Microbial Process for the Production of Benzyl Acetate. Nat Chem Eng 2024. <https://doi.org/10.1038/s44286-023-00022-0>`_

  #. `Characterizing the Opportunity Space for Sustainable Hydrothermal Valorization of Wet Organic Wastes. Environ Sci Technol <https://doi.org/10.1021/acs.est.3c07394>`_

  #. `An End-to-End Pipeline for Succinic Acid Production at an Industrially Relevant Scale Using Issatchenkia Orientalis. Nat Commun 2023 <https://doi.org/10.1038/s41467-023-41616-9>`_

  #. `Metabolic engineering of Yarrowia lipolytica to produce 1,2-diacyl-3-acetyl triacylglycerol. Metabolic Engineering 2023 <https://doi.org/10.1016/j.ymben.2023.01.003>`_

  #. `Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries. ACS Sustainable Chem. Eng. 2022 <https://doi.org/10.1021/acssuschemeng.2c04204>`_

  #. `Adsorptive separation and recovery of triacetic acid lactone from fermentation broth. Biofuels Bioprod Bioref 2022 <https://doi.org/10.1002/bbb.2427>`_

  #. `Rewiring yeast metabolism for producing 2,3-butanediol and two downstream applications: Techno-economic analysis and life cycle assessment of methyl ethyl ketone (MEK) and agricultural biostimulant production. Chemical Engineering Journal 2022 <https://doi.org/10.1016/j.cej.2022.138886>`_

  #. `Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem Eng 2021 <https://doi.org/10.1021/acssuschemeng.1c05441>`_
  
  #. `Renewable linear alpha-olefins by base-catalyzed dehydration of biologically-derived fatty alcohols. Green Chemistry 2021 <https://doi.org/10.1039/D1GC00243K>`_

  #. `Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem Eng 2021 <https://doi.org/10.1021/acssuschemeng.0c08055>`_

  #. `Techno-Economic Evaluation of Biorefineries Based on Low-Value Feedstocks Using the BioSTEAM Software: A Case Study for Animal Bedding. Processes 2020 <https://doi.org/10.3390/pr8080904>`_

* **Plastics and recycling**:

  #. `Screening green solvents for multilayer plastic film recycling processes. Computers & Chemical Engineering 2025 <https://doi.org/10.1016/j.compchemeng.2025.109129>`_

  #. `Comparative Techno-Economic Analysis and Life Cycle Assessment of Producing High-Value Chemicals and Fuels from Waste Plastic via Conventional Pyrolysis and Thermal Oxo-Degradation. Energy Fuels 2023 <https://doi.org/10.1021/acs.energyfuels.3c02321>`_

  #. `High-purity polypropylene from disposable face masks via solvent-targeted recovery and precipitation. Green Chemistry 2023 <https://doi.org/10.1039/D3GC00205E>`_


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`