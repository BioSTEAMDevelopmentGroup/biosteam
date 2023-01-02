# Cane Biorefineries and Benchmarks

The `cane` module contains sugarcane, oilcane, and energycane biorefinery configurations, 
as discussed in [[1]](#1). Two oilcane configurations are currently available: 
(I) oil extraction by the expression of bagasse and centrifugation of vinasse from juice fermentation, and 
(II) oil extraction after an integrated, single-step co-fermentation of both juice 
and bagasse hydrolysate. Mass balances around oil separation are based on experimental results 
and the oil composition in terms of free fatty acids, triacyl glycerides, and polar lipids 
is accounted for. 

```{toctree}
:hidden:
chemicals
systems
units
```

Getting Started
---------------

Four biorefineries can be created using the names detailed in the following table:

| Feedstocks                | Products            | Direct cogeneration | Integrated Co-fermentation|
| ------------------------- | ------------------  | ------------------- | ------------------------- |
| Oilcane                   | Ethanol & biodiesel | O1                  | O2                        |
| Oilcane                   | Ethanol & crude oil | O3                  | O4                        |
| Oilcane & oil-sorghum     | Ethanol & biodiesel | O1\*                | O2\*                      |
| Oilcane & oil-sorghum     | Ethanol & crude oil | O3\*                | O4\*                      |
| Sugarcane                 | Ethanol & biodiesel | S1                  | S2                        |
| Sugarcane & sweet sorghum | Ethanol             | S1\*                | S2\*                      |

Here are a few examples:

```python
>>> from biorefineries import cane
>>> S1 = cane.Biorefinery('S1') # Load conventional sugarcane biorefinery
>>> S1.sys.show(data=False) # Full system
System: sugarcane_sys
Highest convergence error among components in recycle
stream M201-0 after 5 loops:
- flow rate   1.36e-12 kmol/hr (3.8e-14%)
- temperature 4.41e-06 K (1.3e-06%)
ins...
[0] sugarcane
[1] H3PO4
[2] lime
[3] polymer
[4] denaturant
outs...
[0] advanced_ethanol
[1] vinasse
[2] wastewater
[3] emissions
[4] ash_disposal

>>> O1 = cane.Biorefinery('O1') # Load direct cogeneration oilcane configuration
>>> O1.sys.show(data=False)
System: oilcane_sys
Highest convergence error among components in recycle
streams {C504-1, P518-0} after 4 loops:
- flow rate   2.07e-08 kmol/hr (2.1e-06%)
- temperature 6.68e-07 K (2e-07%)
ins...
[0] oilcane
outs...
[0] advanced_ethanol
[1] biodiesel
[2] crude_glycerol
[3] vinasse

```

To retrieve economic and environmental results at different scenarios, you can 
use the Model object:

```python
>>> import biorefineries.cane
>>> import pandas as pd
>>> pd.set_option('display.max_rows', 50) # Show complete data table
>>> O2 = cane.Biorefinery('O2') # Load integrated co-fermentation oilcane configuration
>>> parameters = O2.model.get_baseline_sample() # All parameters at the baseline scenario
>>> parameters
                                   Oil recovery [%]                                 60
                                   Saccharification oil recovery [%]                70
                                   Cane operating days [day/yr]                    180
                                   Sorghum operating days [day/yr]                  45
                                   Annual crushing capacity [MT/yr]            1.6e+06
Stream-cellulosic ethanol          Price [USD/L]                                 0.902
Stream-advanced ethanol            Price [USD/L]                                 0.573
Stream-biodiesel                   Price [USD/L]                                  1.01
Stream-cellulosic based diesel     Price [USD/L]                                   1.5
Stream-natural gas                 Price [USD/m3]                                0.154
                                   Electricity price [USD/kWhr]                 0.0641
                                   IRR [%]                                          10
Stream-crude glycerol              Price [USD/kg]                                 0.16
Stream-pure glycerine              Price [USD/kg]                                 0.65
Saccharification                   Reaction time [hr]                               72
Stream-cellulase                   Price [USD/kg]                                0.212
cellulase                          Cellulase loading [wt. % cellulose]            0.02
Pretreatment reactor system        Base cost [million USD]                    1.97e+07
Pretreatment and saccharification  Cane glucose yield [%]                           91
                                   Sorghum glucose yield [%]                        79
                                   Cane xylose yield [%]                          97.5
                                   Sorghum xylose yield [%]                         86
Cofermenation                      Glucose to ethanol yield [%]                     90
                                   Xylose to ethanol yield [%]                      42
Cofermentation                     Ethanol titer [g/L]                            68.5
                                   Ethanol productivity [g/L]                    0.951
Cofermenation                      Glucose to microbial oil yield [%]             49.5
                                   Xylose to microbial oil yield [%]              49.5
Cofermentation                     Microbial oil titer [g/L]                      27.4
                                   Microbial oil productivity [g/L]               0.31
oilcane                            Cane PL content [% oil]                          10
oilsorghum                         Sorghum PL content [% oil]                       10
oilcane                            Cane FFA content [% oil]                         10
oilsorghum                         Sorghum FFA content [% oil]                      10
oilcane                            Cane oil content [dry wt. %]                     10
oilsorghum                         Relative sorghum oil content [dry wt. %]       -1.5
                                   TAG to FFA conversion [% oil]                    23
Stream-oilcane                     GWP [kg*CO2-eq/kg]                           0.0352
Stream-methanol                    GWP [kg*CO2-eq/kg]                             0.45
Stream-pure glycerine              GWP [kg*CO2-eq/kg]                             1.67
Stream-cellulase                   GWP [kg*CO2-eq/kg]                            0.404
Stream-natural gas                 GWP [kg*CO2-eq/kg]                             0.33
                                   Income tax [%]                                   21
dtype: float64

>>> parameters['oilcane', 'Cane oil content [dry wt. %]'] = 10 # Change oil content
>>> O2.model(parameters) # Evaluate at new oil content
Biorefinery              MFPP [USD/MT]                                          47.1
                         Feedstock consumption [MT/yr]                       1.6e+06
                         Biodiesel production [L/MT]                            27.1
                         Biodiesel yield [L/hc]                                  NaN
                         Ethanol production [L/MT]                              99.8
                         Electricity production [kWhr/MT]                       41.2
                         Net energy production [GGE/MT]                           27
                         Natural gas consumption [m3/MT]                           0
                         TCI [10^6*USD]                                          503
                         Heat exchanger network error [%]                  -1.82e-07
Economic allocation      GWP [kg*CO2*eq / USD]                                 0.483
                         Ethanol GWP [kg*CO2*eq / L]                           0.367
                         Biodiesel GWP [kg*CO2*eq / L]                         0.622
                         Crude glycerol GWP [kg*CO2*eq / kg]                  0.0772
                         Electricity GWP [kg*CO2*eq / MWhr]                     30.9
Displacement allocation  Ethanol GWP [kg*CO2*eq / L]                          -0.446
                         Biodiesel GWP [kg*CO2*eq / L]                         -1.65
Energy allocation        Biofuel GWP [kg*CO2*eq / GGE]                          2.04
                         Ethanol GWP [kg*CO2*eq / L]                            0.36
                         Biodiesel GWP [kg*CO2*eq / L]                         0.566
                         Crude-glycerol GWP [kg*CO2*eq / kg]                   0.216
Biorefinery              IRR [%]                                                8.91
                         MFPP derivative [USD/MT]                               1.33
                         Biodiesel production derivative [L/MT]                 2.71
                         Ethanol production derivative [L/MT]                  -3.05
                         Electricity production derivative [kWhr/MT]            5.94
                         Natural gas consumption derivative [cf/MT]                0
                         TCI derivative [10^6*USD]                              1.93
Economic allocation      GWP derivative [kg*CO2*eq / USD]                   -0.00851
Ethanol                  Ethanol GWP derivative [kg*CO2*eq / L]             -0.00655
Biodiesel                Biodiesel GWP derivative [kg*CO2*eq / L]            -0.0111
Crude glycerol           Crude glycerol GWP derivative [kg*CO2*eq / kg]     -0.00136
Electricity              Electricity GWP derivative [kg*CO2*eq / MWhr]        -0.546
Biorefinery              ROI [%]                                                9.59
                         Competitive oilcane biomass yield [% sugarcane]         NaN
                         Competitive microbial oil yield [wt. %]                 NaN
dtype: float64

```

## References
<a id="1">[1]</a> 
    Cortés-Peña, Y.R., C.V. Kurambhatti, K. Eilts, V. Singh, J.S. Guest, “Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries,” ACS Sustainable Chemistry & Engineering.

