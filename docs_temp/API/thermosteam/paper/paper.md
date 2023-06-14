---
title: "Thermosteam: BioSTEAM's Premier Thermodynamic Engine"
tags:
  - Python
  - BioSTEAM
  - thermodynamics
  - process modeling
  - mass and energy balance
  - chemical properties
  - mixture properties
  - phase equilibrium
authors:
  - name: Yoel Cortés-Peña
    orcid: 0000-0003-1742-5059
    affiliation: "1, 2"
affiliations:
  - name: Department of Civil and Environmental Engineering, University of Illinois at Urbana-Champaign
    index: 1
  - name: DOE Center for Advanced Bioenergy and Bioproducts Innovation (CABBI)
    index: 2
date: 12 May 2020
bibliography: paper.bib
---

# Summary

`Thermosteam` is a thermodynamic engine capable of solving mass and 
energy balances, estimating mixture properties, solving thermodynamic phase 
equilibria, and modeling stoichiometric reactions. All chemical data in 
`Thermosteam` is imported from the `chemicals` library [@chemicals], an 
open-source compilation of data and functions for the estimation of pure 
component chemical and mixture properties. `Thermosteam`'s fast and flexible
platform has enabled the evaluation of conceptual and emerging biochemical production
processes. The Biorefinery Simulation and Techno-Economic Analysis Modules 
(BioSTEAM) — capable of modeling reactors, distillation columns, heat exchangers, 
and other unit operations — has adopted `Thermosteam` as its premier 
thermodynamic engine [@BioSTEAM]. Published biorefinery designs modeled in 
BioSTEAM implement thermodynamic property packages created with `Thermosteam` 
[@Bioindustrial-Park], including a cornstover biorefinery for the production of 
cellulosic ethanol, a lipid-cane biorefinery for the co-production of ethanol 
and biodiesel, and a wheatstraw biorefinery for the production of cellulosic 
ethanol [@BioSTEAM;@Sanchis].

# Statement of Need

The overarching goal of `Thermosteam` is to aid the rigorous design and 
simulation of chemical production processes, whereby low value feedstocks are
converted to high value products via chemical reactions and thermodynamic-driven
separations. For example, modeling the separation of volatile chemicals from 
heavier ones in a distillation column (e.g., distilling ethanol from water),
requires vapor-liquid phase equilibrium calculations to predict how well volatile 
chemicals selectively partition into the vapor phase. Additionally, fluid 
viscosities, densities, and surface tensions are required to appropriately design 
a distillation column that can achieve a specified recovery of chemicals 
[@Perry]. 

Several open-source libraries in Python have comparable capabilities to 
`Thermosteam` in the estimation of fluid properties and phase equilibria: most 
notably `Cantera` and `CoolProp`. `Cantera` is a collection of software tools 
capable of modeling kinetic reactions, thermodynamic equilibrium, and mixture 
properties [@Cantera]. `Cantera`'s built-in chemicals are limited to 8, but new
chemicals can be defined by users with flexibility on the amount of detail. 
`Thermosteam` has yet to implement any features on kinetic reaction networks, 
but exposes a larger set of roughly 20,000 built-in chemicals from the `chemicals` 
library. Users may also define new models and pseudo-chemicals that are compatible 
with all of `Thermosteam`'s features. `CoolProp` offers fast and accurate 
thermodynamic and transport properties for 122 chemical components, and can 
estimate phase equilibrium and mixture properties [@CoolProp]. `CoolProp` 
also offers an interface to the NIST REFPROP software, which is considered the 
gold standard in thermophysical properties [@REFPROP]. It is within 
`Thermosteam`'s roadmap to use `CoolProp` as part of its built-in models. 
While `CoolProp` focuses on thermophysical chemical properties, `Thermosteam` 
also includes mass and energy balances and stoichiometric reactions as one of 
its central features.

# Roadmap

The main development  items in `Thermosteam`'s roadmap concerns the implementation
of fast, robust, and accurate algorithms for estimating mixture properties
and solving thermodynamic phase equilibria. Through `Thermosteam`, `BioSTEAM` is
able to evaluate a range of biofuels and bioproducts, but further efforts on these 
development items would enable the evaluation of a broader portfolio of 
potential bioproducts. 

In `Thermosteam`, Peng Robinson is the default equation of state 
for all pure components. However, the estimation of pure component chemical 
properties is not limited to solving the equation of state. Several models 
of thermodynamic properties (e.g., density, heat capacity, vapor pressure, 
heat of vaporization) are correlations that rely on fitted coefficients 
and key chemical properties (e.g., critical temperature and pressure). To 
facilitate the calculation of mixture properties, `Thermosteam`'s 
mixing rule estimates mixture properties by assuming a molar weighted average 
of the pure chemical properties. However, `Thermosteam` aims to implement rigorous 
equation of state (EOS) mixing rules for the estimation of mixture properties.

`Thermosteam` allows for fast estimation of thermodynamic equilibrium within 
hundreds of microseconds through the smart use of cache and Numba just-in-time 
(JIT) compiled functions [@numba]. The main vapor-liquid equilibrium (VLE) 
algorithm solves the modified Raoult’s law equation with activity coefficients
estimated through UNIQUAC Functional-group Activity Coefficients (UNIFAC) 
interaction parameters [@Gmehling]. Modified Raoult’s law is suitable to 
estimate VLE of nonideal mixtures under low to moderate pressures. At high to 
near-critical pressures, gaseous nonidealities become more significant. In a 
near future, `Thermosteam` may also implement the Predictive Soave–Redlich–Kwong
(PSRK) functional group model for estimating phase equilibrium of critical
mixtures. 

All of `Thermosteam`'s application program interface (API) is documented with 
examples. These examples also serve as preliminary tests that must pass before
accepting any changes to the software via continuous integration on Github.
Additionally, the online documentation includes a full tutorial that concludes 
with the creation of a property package. `Thermosteam`’s powerful features 
and extensive documentation encourage users to become a part of its
community-driven platform and help it become more industrially and academically 
relevant. 

# Acknowledgements

I would like to thank Caleb Bell for developing the open-source `chemicals` library
in Python, which has served as both groundwork and inspiration for developing `Thermosteam`. 
This material is based upon work supported by the National Science Foundation Graduate
Research Fellowship Program under Grant No. DGE—1746047. Any opinions, findings,
and conclusions or recommendations expressed in this publication are those of 
the authors and do not necessarily reflect the views of the National Science
Foundation. This work was funded by the DOE Center for Advanced Bioenergy and 
Bioproducts Innovation (U.S. Department of Energy, Office of Science, Office of 
Biological and Environmental Research under Award Number DE-SC0018420). Any 
opinions, findings, and conclusions or recommendations expressed in this 
publication are those of the author and do not necessarily reflect the views of
the U.S. Department of Energy.

# References