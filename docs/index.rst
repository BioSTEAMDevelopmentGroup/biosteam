Fast and Flexible Process Simulation
====================================

.. toctree::
   :maxdepth: 2
   :hidden:
   
   tutorial/index
   API/index
   impact/index
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
        analysis (TEA) and life-cycle assessment (LCA) of chemical processes across thousands 
        of scenarios. It features over 20,000 chemicals in the databank, 
        rigorous thermodynamic property packages, and advanced unit operations. 
        Users can also define custom chemicals, unit operations, and process specifications.
        Given the open nature of BioSTEAM, it serves as a platform to develop accessible
        process models, build specialized simulation tools, 
        and advance process simulation algorithms.
        

.. grid:: 1 2 3 4
  
    .. grid-item-card:: Getting Started
       :text-align: center
       :link: https://biosteam.readthedocs.io/en/latest/tutorial/Getting_started.html
       :link-type: url
       :padding: 1

       .. image:: images/extraction.png
          :height: 100
          :class: only-light
          :align: center
          
       .. image:: images/extraction_dark.png
          :height: 100
          :class: only-dark
          :align: center
          
       Comprehensive tutorials


    .. grid-item-card:: User Manual
       :text-align: center
       :link: https://biosteam.readthedocs.io/en/latest/API/index.html
       :link-type: url
       :padding: 1
       
       .. image:: images/API_logo.png
          :height: 100
          :class: only-light
          :align: center
          
       .. image:: images/API_logo_dark.png
          :height: 100
          :class: only-dark
          :align: center
       
       Detailed documentation

       
    .. grid-item-card:: Bioindustrial-Park
       :text-align: center
       :link: https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park
       :link-type: url
       :padding: 1
          
       .. image:: images/membrane_bioreactor.png
          :height: 100
          :class: dark-light
          :align: center
          
       Process models
       
       
    .. grid-item-card:: YouTube Channel
       :text-align: center
       :link: https://www.youtube.com/@yoelcortes-pena2100/videos
       :link-type: url
       :padding: 1
    
       .. image:: images/analysis.png
          :height: 100
          :class: only-light
          :align: center
          
       .. image:: images/analysis_dark.png
          :height: 100
          :class: only-dark
          :align: center
    
       Recorded workshops


Installation
------------

#. If you have an installation of Python with pip, simple install it with:

   .. code-block:: bash

      $ pip install biosteam


   To get the git version, run:

   .. code-block:: bash
   
      $ git clone --depth 10 git://github.com/BioSTEAMDevelopmentGroup/biosteam
      
   If you are a developer, we recommend getting BioSTEAM together with the BioIndustrial-Park and ThermoSTEAM:
   
      $ git clone --recurse-submodules -j8 --depth 10 git://github.com/BioSTEAMDevelopmentGroup/biosteam

#. BioSTEAM uses `Graphviz <http://www.graphviz.org/>`__ to make flowsheet diagrams. To properly install Graphviz in an anaconda distribution, run the following line:

   .. code-block:: bash
    
      $ conda install python-graphviz

#. If you **do not** have an anaconda distribution, you will need to install Graphviz separately as follows:

   * Windows: Download the EXE installer and follow the instructions listed `in this link <https://graphviz.org/download/>`__

   * Ubuntu: 

     .. code-block:: bash
    
       $ sudo apt-get install graphviz
   
   * MacOS: 

     .. code-block:: bash
    
        $ brew install graphviz


Common Issues
-------------

* **Unit and system diagrams are not displaying:**

  Graphviz may not be properly installed or may be missing from your python path. 
  Please follow the graphviz installation procedure outlined above.

* **Cannot install/update BioSTEAM:**

  If you are having trouble installing or updating BioSTEAM, it may be due to dependency issues. 
  Some dependencies like chaospy/numpoly require Microsoft C++ Build Tools. Download and install the `C++ build tools here. <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`__
  You can also bypass dependency issues using:
  
  .. code-block:: bash

     $ pip install --user --ignore-installed biosteam

  You can make sure you install the right version by including the version number:

  .. code-block:: bash

     $ pip install biosteam==<version>
