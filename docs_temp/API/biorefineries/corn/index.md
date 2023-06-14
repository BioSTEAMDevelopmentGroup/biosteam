# Corn Conventional dry grind process

The `corn` module includes systems and units for the co-production of ethanol, distiller's dried grains with solubles (DDGS), and low-grade oil from corn kernels through the conventional dry grind process,
[[1]](#1),[[2]](#2) the biorefinery is modeled based on Kwiatkowski et al. [[1]](#1).

```{toctree}
:hidden:
chemicals
systems
units
```

Getting Started
---------------
```python
>>> from biorefineries import corn
>>> cn = corn.Biorefinery('conventional dry-grind')
>>> cn.corn_sys.show()
System: corn_sys
Highest convergence error among components in recycle
stream S1-0 after 4 loops:
- flow rate   9.61e+00 kmol/hr (0.16%)
- temperature 7.77e-02 K (0.021%)
ins...
[0] alpha_amylase
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water           1.53
                    SolubleProtein  0.000139
[1] ammonia
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): NH3  4.61
[2] lime
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): CaO  0.0841
[3] recycled_process_water
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water  3.88e+03
[4] steam
    phase: 'g', T: 508.99 K, P: 3.11e+06 Pa
    flow (kmol/hr): Water  911
[5] gluco_amylase
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water           4.36
                    SolubleProtein  0.000533
[6] sulfuric_acid
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): H2SO4  0.4
[7] yeast
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water  3.52
                    Yeast  0.123
[8] scrubber_water
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water  1.01e+03
[9] denaturant
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Octane  5.73
[10] makeup_water
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow: 0
[11] s83
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow: 0
[12] dryer_air
    phase: 'g', T: 298.15 K, P: 1.01325e+06 Pa
    flow (kmol/hr): O2  285
                    N2  794
[13] natural_gas
    phase: 'g', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): CH4  66.8
[14] oxidizer_air
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow: 0
[15] s82
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): CH4  3.93
[16] corn
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water             385
                    Ash               647
                    TriOlein          1.77
                    Starch            174
                    Fiber             30.4
                    SolubleProtein    9.69
                    InsolubleProtein  14.1
[17] s84
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow: 0
outs...
[0] s41
    phase: 'g', T: 303.85 K, P: 101325 Pa
    flow (kmol/hr): CO2  330
[1] ethanol
    phase: 'l', T: 338.61 K, P: 101325 Pa
    flow (kmol/hr): Water    6.9
                    Ethanol  324
                    Octane   5.73
[2] crude_oil
    phase: 'l', T: 342.97 K, P: 30953 Pa
    flow (kmol/hr): TriOlein  0.774
[3] process_water
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water  3.88e+03
[4] wastewater
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water  946
[5] s79
    phase: 'g', T: 373.15 K, P: 101325 Pa
    flow (kmol/hr): Water  134
                    CO2    66.8
[6] s80
    phase: 'l', T: 637.19 K, P: 101325 Pa
    flow (kmol/hr): Water  1.14e+03
                    CO2    3.94
                    O2     277
                    N2     794
[7] s4
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (kmol/hr): Water             1.15
                    Ash               1.94
                    TriOlein          0.00532
                    Starch            0.523
                    Fiber             0.0912
                    SolubleProtein    0.0291
                    InsolubleProtein  0.0422
[8] DDGS
    phase: 'l', T: 343.15 K, P: 101325 Pa
    flow (kmol/hr): Water             65.5
                    Ash               645
                    Yeast             12
                    CaO               0.0841
                    TriOlein          0.995
                    H2SO4             0.4
                    Fiber             30.3
                    SolubleProtein    9.66
                    InsolubleProtein  14
[9] s81
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow: 0
```

References
----------
<a id="1">[1]</a> 
    J.R. Kwiatkowski, A. J. McAloon, F. Taylor and D. B. Johnston,
    Modeling the process and costs of fuel ethanol production by the corn 
    dry-grind process, Industrial Crops and Products, 2006, 23, 288–296. 
    https://www.sciencedirect.com/science/article/pii/S0926669005000944

<a id="2">[2]</a> 
    C. Kurambhatti, D. Kumar and V. Singh, Impact of Fractionation Process 
    on the Technical and Economic Viability of Corn Dry Grind Ethanol Process,
    Processes, 2019, 7, 578. https://www.mdpi.com/2227-9717/7/9/578