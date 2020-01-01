# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 02:13:07 2019

@author: yoelr
"""
from .utils import assert_unit_streams, assert_unit_results

__all__ = ('test_binary_distillation',)

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
        '                   Flow                      kmol/hr  3.33e+03\n'
        '                   Cost                       USD/hr      1.62\n'
        'Low pressure steam Duty                        kJ/hr  9.03e+06\n'
        '                   Flow                      kmol/hr       232\n'
        '                   Cost                       USD/hr      55.2\n'
        'Design             Theoretical feed stage                    9\n'
        '                   Theoretical stages                       13\n'
        '                   Minimum reflux              Ratio     0.687\n'
        '                   Reflux                      Ratio      1.37\n'
        '                   Rectifier stages                         14\n'
        '                   Stripper stages                          13\n'
        '                   Rectifier height               ft      33.2\n'
        '                   Stripper height                ft      31.7\n'
        '                   Rectifier diameter             ft      7.83\n'
        '                   Stripper diameter              ft         3\n'
        '                   Rectifier wall thickness       in     0.312\n'
        '                   Stripper wall thickness        in      0.25\n'
        '                   Rectifier weight               lb  9.94e+03\n'
        '                   Stripper weight                lb  3.31e+03\n'
        'Cost               Rectifier trays               USD  2.58e+04\n'
        '                   Stripper trays                USD  1.21e+04\n'
        '                   Rectifier tower               USD  1.16e+05\n'
        '                   Stripper tower                USD  6.09e+04\n'
        '                   Condenser                     USD  3.34e+04\n'
        '                   Boiler                        USD   2.6e+04\n'
        'Purchase cost                                    USD  2.74e+05\n'
        'Utility cost                                  USD/hr      56.9')