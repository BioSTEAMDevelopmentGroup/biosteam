# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import pytest

def test_dryer_relative_humidity_and_energy_balance():
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
    p_sat_water = hot_air.thermo.chemicals.Water.Psat(hot_air.T)
    rh = p_water/p_sat_water

    assert 0 <= rh <= 1, f"Relative humidity must be between 0 and 1. Current: {rh:.2f}"
    assert p_sat_water > 0, "Psat must be > 0."
    assert rh == pytest.approx(D1.RH, rel=1e-3), f"Relative humidity of hot air ({rh}) did not match DrumDryer RH ({D1.RH})."

    # Energy balance
    ## Water evaporated
    dried_solids = D1.outs[0]

    water_in = feed.imass['Water']
    water_solids = dried_solids.imass['Water']
    water_evaporated = water_in - water_solids

    water_hot_air = hot_air.imass['Water']

    assert water_evaporated > 0, "No water was evaporated"
    assert water_evaporated == pytest.approx(hot_air.imass['Water']), f"water_evaporated ({water_evaporated:.2f} kg/hr) do not match water in hot_air ({water_hot_air:.2f} kg/hr)"

    ## Hvap in kJ/kg
    Hvap_mol = D1.thermo.chemicals.Water.Hvap(hot_air.T)
    MW_water = D1.thermo.chemicals.Water.MW
    Hvap_kJ_per_kg = Hvap_mol / MW_water

    ## Water evaporation duty
    Q_evap_required = water_evaporated * Hvap_kJ_per_kg
    
    assert Q_evap_required > 0, "Evaporation duty must be positive"

    ## Methane duty
    Q_methane = D1.ins[2].F_mol * D1.chemicals.CH4.LHV

    assert Q_methane > Q_evap_required, (
        f"Methane energy ({Q_methane:.2f} kJ/hr) is lower than"
        f"theoretical evaporation duty ({Q_evap_required:.2f} kJ/hr)"
    )
    
    ## Energy balance
    H_in = sum(s.H for s in D1.ins)
    H_out = sum(s.H for s in D1.outs)
    Q_enthalpy = H_out - H_in
    
    assert Q_methane == pytest.approx(Q_enthalpy, rel=0.005), (
        f"Global energy balance is not consistent:\n"
        f"Q_enthalpy = {Q_enthalpy:.2f} kJ/hr\n"
        f"Q_methane = {Q_methane:.2f} kJ/hr\n"
        f"ratio enthalpy/methane = {Q_enthalpy / Q_methane:.4f}"
    )

if __name__ == '__main__':
    test_dryer_relative_humidity_and_energy_balance()