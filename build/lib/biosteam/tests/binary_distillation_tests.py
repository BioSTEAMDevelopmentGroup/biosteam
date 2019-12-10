# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:13:07 2019

@author: yoelr
"""
from .utils import assert_unit_streams, assert_unit_results

def test_binary_distillation():
    from biosteam import Species, Stream
    from biosteam.units import Distillation
    
    # Set up stream
    Stream.species = Species('Water', 'Methanol', 'Glycerol')
    feed = Stream('feed', flow=(80, 100, 25))
    feed.T = feed.bubble_T()[0] # Feed at bubble point T
    
    # Set up column
    D1 = Distillation('D1', ins=feed, outs=('distillate', 'bottoms'),
                      LHK=('Methanol', 'Water'),
                      y_top=0.99, x_bot=0.01, k=2)
    D1.is_divided = True
    D1.simulate()
    
    assert_unit_streams(D1, 
        "Distillation: D1\n"
        "ins...\n"
        "[0] feed\n"
        "    phase: 'l', T: 349.28 K, P: 101325 Pa\n"
        "    flow (kmol/hr): Water     80\n"
        "                    Methanol  100\n"
        "                    Glycerol  25\n"
        "outs...\n"
        "[0] distillate\n"
        "    phase: 'g', T: 338.06 K, P: 101325 Pa\n"
        "    flow (kmol/hr): Water     1\n"
        "                    Methanol  99.2\n"
        "[1] bottoms\n"
        "    phase: 'l', T: 373.21 K, P: 101325 Pa\n"
        "    flow (kmol/hr): Water     79\n"
        "                    Methanol  0.798\n"
        "                    Glycerol  25")
    assert_unit_results(
        'Distillation                                   Units        D1\n'
        'Cooling water      Duty                        kJ/hr -4.86e+06\n'
        '                   Flow                      kmol/hr  7.18e+03\n'
        '                   Cost                       USD/hr         0\n'
        'Low pressure steam Duty                        kJ/hr  1.01e+07\n'
        '                   Flow                      kmol/hr       334\n'
        '                   Cost                       USD/hr       102\n'
        'Design             Theoretical feed stage                    4\n'
        '                   Theoretical stages                       13\n'
        '                   Minimum reflux              Ratio     0.687\n'
        '                   Reflux                      Ratio      1.37\n'
        '                   Rectifier stages                          6\n'
        '                   Stripper stages                          27\n'
        '                   Rectifier height               ft      21.4\n'
        '                   Stripper height                ft      52.4\n'
        '                   Rectifier diameter             ft      7.82\n'
        '                   Stripper diameter              ft         3\n'
        '                   Rectifier wall thickness       in     0.312\n'
        '                   Stripper wall thickness        in      0.25\n'
        '                   Rectifier weight               lb  6.96e+03\n'
        '                   Stripper weight                lb  5.31e+03\n'
        'Cost               Rectifier trays               USD  1.26e+04\n'
        '                   Stripper trays                USD  1.75e+04\n'
        '                   Rectifier tower               USD  8.07e+04\n'
        '                   Stripper tower                USD  7.81e+04\n'
        '                   Condenser                     USD  2.62e+04\n'
        '                   Boiler                        USD   2.2e+04\n'
        'Purchase cost                                    USD  2.37e+05\n'
        'Utility cost                                  USD/hr       102')