===============================================================================
BioSTEAM: Los Modulos de Simulación y Análisis Tecno-Económico de Biorefinerias
===============================================================================

.. image:: http://img.shields.io/pypi/v/biosteam.svg?style=flat
   :target: https://pypi.python.org/pypi/biosteam
   :alt: Versión
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://biosteam.readthedocs.io/en/latest/
   :alt: Documentación
.. image:: http://img.shields.io/badge/license-UIUC-blue.svg?style=flat
   :target: https://github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
   :alt: Licencia
.. image:: https://img.shields.io/pypi/pyversions/biosteam.svg
   :target: https://pypi.python.org/pypi/biosteam
   :alt: Versiones_compatibles
.. image:: https://zenodo.org/badge/164639830.svg
   :target: https://zenodo.org/badge/latestdoi/164639830
.. image:: https://coveralls.io/repos/github/BioSTEAMDevelopmentGroup/biosteam/badge.svg?branch=master
   :target: https://coveralls.io/github/BioSTEAMDevelopmentGroup/biosteam?branch=master
.. image:: https://badges.gitter.im/BioSTEAM-users/BioSTEAM.svg
   :alt: Únete al chat https://gitter.im/BioSTEAM-users/community
   :target: https://gitter.im/BioSTEAM-users/community

**Lea en:** `English <README.rst>`_

.. contents::

¿Qué es BioSTEAM?
-----------------

BioSTEAM es un projecto para llevar a cabo el diseño, simulación, 
análisis tecno-económico, y análisis de cyclo de vida de biorrefinerías bajo
incertidumbre [1]_. BioSTEAM esta diseñada para facilitar la toma de decisiones
al simular múltiples escenarios que podrían surgir dentro del desarrollo e
implementación de una biorrefinería. Configuraciones de biorrefinerías
estan disponibles en el `Bioindustrial-Park <https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park>`_, 
el reposito principal de BioSTEAM para modelos y resultados. El crecimiento y 
mantenimiento de BioSTEAM a largo plazo esta apollada por ambos la communidad de 
dessarrollo de codigo abierto y las instituciones auspiciando a BioSTEAM, incluyendo
el `Centro para Innovaciones de Bioenergia y Bioproductos Avanzados (CABBI) <https://cabbi.bio/>`_. 
Mediante su plataforma de dessarrollo abierto, BioSTEAM aspira fomentar
communicacion y transparencia dentro de la communidad de dessarrollo de biorefinerias 
para un esfuerzo integrado en acelerar la evaluación de biocombustibles y 
bioproductos.

Toda data en quimicos y algoritmos para estimar propiedades quimicos son importados 
de `chemicals <https://github.com/CalebBell/chemicals>`_ y `thermo <https://github.com/CalebBell/chemicals>`_,
paquetes abiertos dessarollado por Caleb Bell. El motor thermodynamico principal de BioSTEAM, 
`ThermoSTEAM <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_, facilita
la creación de paquetes de propiedades thermodynamicas atravéz de estos paquetes.

Instalación
-----------

Encuentra la ultima versión de BioSTEAM en `PyPI <https://pypi.python.org/pypi/biosteam/>`_. 
Si tienes una instalación de Python con pip, instalalo con:

    $ pip install biosteam

Para la versión git, corra:

    $ git clone git://github.com/BioSTEAMDevelopmentGroup/biosteam

Para ayuda con la instalación, visita la `documentación <https://biosteam.readthedocs.io/en/latest/#installation>`_.

Documentación
-------------

La documentación de BioSTEAM esta disponible en la web (por ahora, solo en inglés):

    http://biosteam.readthedocs.io/

Reporte de Errores 
------------------

Para reportar errores, por favor publique una entrada en:

    https://github.com/BioSTEAMDevelopmentGroup/biosteam

Contribuciones
--------------

Para instrucciones sobre cómo contribuir, visite:

    https://biosteam.readthedocs.io/en/latest/contributing/index.html


Licencia
--------

Vea ``LICENSE.txt`` para información sobre los terminos & condiciones del uso
de esta aplicación, y una EXLUSIÓN DE TODA GUARANTÍA.

Aunque no sea requerido por la licencia de BioSTEAM, si te conviene,
por favor cite a BioSTEAM si lo utilizaste. Para el beneficio de la communidad, 
tambien considere contribuyir cualquier cambio o mejora que hallas hecho.


Sobre los autores
-----------------

BioSTEAM fue creada y dessarrollada por `Yoel Cortés-Peña <https://yoelcortes.github.io/me/>`_ 
como parte de el `Grupo Guest <http://engineeringforsustainability.com/yoelcortespena>`_ 
y el `Centro para Innovaciones de Bioenergia y Bioproductos Avanzados (CABBI) <https://cabbi.bio/>`_ 
en la `Universidad de Illinois en Urbana-Champaign (UIUC) <https://illinois.edu/>`_. 

Referencias
-----------
.. [1] `Cortés-Peña et al. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020. <https://doi.org/10.1021/acssuschemeng.9b07040>`_.


