# cane: Cane Biorefineries and Benchmarks

This module contains sugarcane, oilcane, and energycane biorefinery configurations, 
as discussed in [[1]](#1). Two oilcane configurations are currently available: 
(I) oil extraction by the expression of bagasse and centrifugation of vinasse from juice fermentation, and 
(II) oil extraction after an integrated, single-step co-fermentation of both juice 
and bagasse hydrolysate. Mass balances around oil separation are based on experimental results 
and the oil composition in terms of free fatty acids, triacyl glycerides, and polar lipids 
is accounted for. 

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
>>> import biorefineries.cane
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

>>> oc.load('O1') # Load direct cogeneration oilcane configuration
>>> oc.sys.show(data=False)
System: oilcane_sys
Highest convergence error among components in recycle
streams {C604-1, P618-0} after 4 loops:
- flow rate   2.07e-08 kmol/hr (2.1e-06%)
- temperature 1.00e-06 K (3.1e-07%)
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
>>> O2 = cane.Biorefinery('O2') # Load integrated co-fermentation oilcane configuration
>>> parameters = O2.model.get_baseline_sample() # All parameters at the baseline scenario
>>> parameters
biorefinery                        Oil recovery [%]                                 60
                                   Saccharification oil recovery [%]                70
                                   Cane operating days [day/yr]                    180
                                   Sorghum operating days [day/yr]                  45
                                   Annual crushing capacity [MT/yr]            1.6e+06
Stream-cellulosic ethanol          Price [USD/L]                                 0.902
Stream-advanced ethanol            Price [USD/L]                                 0.573
Stream-biodiesel                   Price [USD/L]                                  1.01
Stream-natural gas                 Price [USD/m3]                                0.123
biorefinery                        Electricity price [USD/kWhr]                 0.0641
                                   IRR [%]                                          10
Stream-crude glycerol              Price [USD/kg]                                 0.16
Stream-pure glycerine              Price [USD/kg]                                 0.65
Saccharification                   Reaction time [hr]                               72
cellulase                          Price [USD/kg]                                0.212
                                   Cellulase loading [wt. % cellulose]            0.02
Pretreatment reactor system        Base cost [million USD]                    1.97e+07
Pretreatment and saccharification  Cane glucose yield [%]                           91
                                   Sorghum glucose yield [%]                        79
                                   Cane xylose yield [%]                          97.5
                                   Sorghum xylose yield [%]                         86
Cofermenation                      Glucose to ethanol yield [%]                     90
                                   Xylose to ethanol yield [%]                      42
Cofermentation                     Titer [g/L]                                    68.5
                                   Productivity [g/L]                            0.951
oilcane                            Cane PL content [% oil]                          10
oilsorghum                         Sorghum PL content [% oil]                       10
oilcane                            Cane FFA content [% oil]                         10
oilsorghum                         Sorghum FFA content [% oil]                      10
oilcane                            Cane oil content [dry wt. %]                     10
oilsorghum                         Relative sorghum oil content [dry wt. %]       -1.5
biorefinery                        TAG to FFA conversion [% oil]                    23
Stream-oilcane                     GWP [kg*CO2-eq/kg]                           0.0352
Stream-methanol                    GWP [kg*CO2-eq/kg]                             0.45
Stream-pure glycerine              GWP [kg*CO2-eq/kg]                             1.67
Stream-cellulase                   GWP [kg*CO2-eq/kg]                            0.161
Stream-natural gas                 GWP [kg*CO2-eq/kg]                             0.33
biorefinery                        Income tax [%]                                   21
dtype: float64

>>> parameters['oilcane', 'Cane oil content [dry wt. %]'] = 10 # Change oil content
>>> oc.model(parameters) # Evaluate at new oil content
Biorefinery              MFPP [USD/MT]                                         42.5
                         Feedstock consumption [MT/yr]                      1.6e+06
                         Biodiesel production [L/MT]                             26
                         Ethanol production [L/MT]                             93.1
                         Electricity production [kWhr/MT]                     0.203
                         Natural gas consumption [m3/MT]                          0
                         TCI [10^6*USD]                                         466
                         Heat exchanger network error [%]                  2.12e-07
Economic allocation      GWP [kg*CO2*eq / USD]                                0.452
                         Ethanol GWP [kg*CO2*eq / L]                          0.339
                         Biodiesel GWP [kg*CO2*eq / L]                        0.456
                         Crude glycerol GWP [kg*CO2*eq / kg]                 0.0723
                         Electricity GWP [kg*CO2*eq / MWhr]                      29
Displacement allocation  Ethanol GWP [kg*CO2*eq / L]                         -0.407
Energy allocation        Biofuel GWP [kg*CO2*eq / GGE]                         1.81
                         Ethanol GWP [kg*CO2*eq / L]                          0.319
                         Biodiesel GWP [kg*CO2*eq / L]                        0.501
                         Crude-glycerol GWP [kg*CO2*eq / kg]                  0.192
Biorefinery              MFPP derivative [USD/MT]                              1.28
                         Biodiesel production derivative [L/MT]                 2.6
                         Ethanol production derivative [L/MT]                 -3.29
                         Electricity production derivative [kWhr/MT]           10.8
                         Natural gas consumption derivative [cf/MT]               0
                         TCI derivative [10^6*USD]                             3.16
Economic allocation      GWP derivative [kg*CO2*eq / USD]                  -0.00407
Ethanol                  Ethanol GWP derivative [kg*CO2*eq / L]            -0.00309
Biodiesel                Biodiesel GWP derivative [kg*CO2*eq / L]          -0.00411
Crude glycerol           Crude glycerol GWP derivative [kg*CO2*eq / kg]   -0.000651
Electricity              Electricity GWP derivative [kg*CO2*eq / MWhr]       -0.261
dtype: float64

```

## References
<a id="1">[1]</a> 
    Cortés-Peña, Y.R., C.V. Kurambhatti, K. Eilts, V. Singh, J.S. Guest, “Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries,” ACS Sustainable Chemistry & Engineering.

