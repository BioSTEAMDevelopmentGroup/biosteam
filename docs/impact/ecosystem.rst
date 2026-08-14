Ecosystem
=========

The BioSTEAM platform
---------------------

Data on chemicals and algorithms to estimate thermodynamic properties are 
imported from `chemicals <https://github.com/CalebBell/chemicals>`__
and `thermo <https://github.com/CalebBell/chemicals>`__,
community-driven open-source libraries developed by Caleb Bell. BioSTEAM's 
premire thermodynamic engine, `ThermoSTEAM <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`__, 
builds upon these libraries to facilitate the creation of thermodynamic property packages.

The `BioSTEAM repository <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`__ houses the unit operations, system convergence algorithms,
the techno-economic/lifecycle assessment frameworks, and the uncertainty/sensitivity analysis features.
The unit operation models use standard modeling algorithms, design procedures, and cost correlations adapted 
from textbooks and literature.

Detailed process models can be found in the `Bioindustrial-Park <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`__ GitHub repository, 
BioSTEAM's premier repository for models and results. You can also find there specialized 
unit operation models and tailored TEA/LCA frameworks.
While the standard for any contribution to BioSTEAM and ThermoSTEAM are high 
(demanding extensive tests and documentation for each new feature), the requirement for
contributing a process model to the BioIndustrial-Park is low to encourage users
to upload models featured in scientific literature.

External tools
--------------

Several groups are leveraging BioSTEAM (and open-source models developed in BioSTEAM)
to build new tools, some of which have AI-integration. Tools integrated with AI
are still in early stages of development. Here are a few of the projects which 
have been openly shared.


.. card:: 
   :link: https://qsdsan.com/
   :link-alt: qsdsan.com

   .. image:: ../images/qsdsan.png
      :height: 100
      :class: only-light
      :align: center
      
   .. image:: ../images/qsdsan_dark.png
      :height: 100
      :class: only-dark
      :align: center
      
   QSDsan is a platform for quantitative sustainable design of sanitation and 
   resource recovery systems led by Prof. Yalin Li at Rutgers University.
   It leverages BioSTEAM's core simulation capabilities.


.. card::
   :link: https://github.com/sarangbhagwat/nskinetics
   :link-alt: github.com/sarangbhagwat/nskinetics
   
   .. image:: ../images/NSKinetics.png
      :height: 100
      :class: only-light
      :align: center
      
   .. image:: ../images/NSKinetics_dark.png
      :height: 100
      :class: only-dark
      :align: center
      
   NSKinetics is a fast, flexible, and convenient package in Python for 
   simulating steady- and non-steady-state reaction kinetics which
   can then drive a BioSTEAM unit operations. It is led by Dr. Sarang Bhagwat.


.. card::
   :link: https://projectpisces.org/
   :link-alt: projectpisces.org
   
   .. image:: ../images/pisces_icon.png
      :height: 100
      :class: only-light
      :align: center
      
   .. image:: ../images/pisces_icon_dark.png
      :height: 100
      :class: only-dark
      :align: center
      
   PISCES harmonizes process data from diverse sources into a machine-readable 
   format to enable AI. It is led by Corrine Scown as part of the Lawrence Berkely National Lab.
   It leverages process models from the BioIndustrial-Park.
   

.. card:: 
   :link: https://puranwater.com/
   :link-alt: puranwater.com
   
   .. image:: ../images/Puran_icon.png
      :height: 100
      :class: only-light
      :align: center
      
   .. image:: ../images/Puran_icon_dark.png
      :height: 100
      :class: only-dark
      :align: center
      
   Puran Water designs industrial wastewater treatment systems and supplies the 
   critical equipment that anchors them. It is led by Hersh Kshetry as part of Puran Water LLC.
   It leverages models built in QSDsan.
   

.. card::
   :link: https://github.com/youmustfight/ai-teas
   :link-alt: github.com/youmustfight/ai-teas
   
   .. image:: ../images/homeworld_collective.png
      :height: 100
      :class: only-light
      :align: center
      
   .. image:: ../images/homeworld_collective_dark.png
      :height: 100
      :class: only-dark
      :align: center
   
   This AI-TEA tool uses LLMs parse existing technoeconomic analysis literature to 
   create process models in BioSTEAM. It is led by Homeworld Collective,
   a US non-profit.
   
   
      
   