# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from numpy.testing import assert_allclose

def test_dryer_relative_humidity():
    settings = bst.settings
    settings.set_thermo(['Water','Protein','O2','N2','CO2','CH4'], cache=True)
    settings.mixture.include_excess_energies = True

    feed = bst.Stream(None,Water=1000,Protein=200,units='kg/hr')
    D1 = bst.DrumDryer(None,(feed,),split={"Protein":0.0})
    D1.simulate()

    # Determine relative humidity of hot air
    hot_air = D1.outs[1]
    
    ## Water partial pressure: p_water = (n_water/n_total) * P 
    n_water = hot_air.imol['Water']
    n_total = hot_air.F_mol
    p_total = hot_air.P
    
    assert n_total > 0, "Hot gas stream is empty"
    
    y_water = (n_water/n_total)
    p_water = y_water * p_total

    ## Relative humidity: rh = p_water / p_saturation
    p_sat_water = settings.chemicals.Water.Psat(hot_air.T)
    rh = p_water/p_sat_water

    assert 0 >= rh >= 1, f"Relative humidity must be between 0 and 1. Current: {rh:.2f}"
    assert p_sat_water > 0, "Psat must be > 0."

if __name__ == '__main__':
    test_dryer_relative_humidity()