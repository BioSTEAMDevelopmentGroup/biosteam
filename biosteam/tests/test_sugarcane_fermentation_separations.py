# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 00:16:24 2020

@author: yrc2
"""
from biosteam import units
import biosteam as bst
from thermosteam import functional as fn
import thermosteam as tmo

__all__ = ('test_sugarcane_fermentation_separations',
)

def test_sugarcane_fermentation_separations():
    return 
    ### Create property package ###
    
    chemicals = tmo.Chemicals(
        ['Water', 'Ethanol', 'Glucose',
         'Sucrose', 'H3PO4', 'P4O10',
         'CO2', 'Octane', 'O2']
    )
    
    Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2 = chemicals
    
    CO2.at_state(phase='g')
    H3PO4.at_state(phase='s')
    P4O10.at_state(phase='s')
    Glucose.at_state(phase='s')
    Sucrose.at_state(phase='s')
    
    DryYeast = tmo.Chemical('DryYeast', MW=1., phase='s',
                            search_db=False, default=True, CAS='Yeast')
    chemicals.append(DryYeast)
    
    Ash = tmo.Chemical('Ash', MW=1., search_db=False, phase='s', default=True)
    chemicals.append(Ash)
    
    # Insolubles occupy a significant volume
    insoluble_solids = (Ash, DryYeast)
    
    # Solubles don't occupy much volume
    soluble_solids = (H3PO4, Glucose, Sucrose) 
    
    for chemical in insoluble_solids:
        V = fn.rho_to_V(rho=1540, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    for chemical in soluble_solids:
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    # Add constant models for molar heat capacity of solids
    Ash.Cn.add_model(0.09 * 4.184 * Ash.MW) 
    
    for chemical in chemicals: chemical.default()
    bst.settings.set_thermo(chemicals)
    chemicals.set_synonym('Water', 'H2O')
    
    ### Create fresh streams ###
    
    # Fresh water
    stripping_water = bst.Stream('stripping_water', Water=5000, units='kg/hr')
    sugar_solution = bst.Stream('sugar_solution', phase='l',
                                Glucose=21.11,
                                Sucrose=125.9,
                                Water=1370 + 4409,
                                H3PO4=0.7575,
                                DryYeast=1.03e+04)
    
    # Ethanol Production
    R301 = units.Fermentation('R301', outs=('CO2', ''), tau=9, efficiency=0.90, N=4) 
    T301 = units.StorageTank('T301', tau=4, vessel_material='Carbon steel')
    T301.line = 'Beer tank'
    
    D301 = units.VentScrubber('D301', ins=(stripping_water, R301-0), gas=('CO2',))
    
    # Separate 99% of yeast
    C301 = units.SolidsCentrifuge('C301', outs=('', 'recycle_yeast'),
                                split=(1, 0.99999, 1, 0.96, 0.01),
                                order=('Ethanol', 'Glucose', 'H3PO4', 
                                       'Water', 'DryYeast'),
                                solids=('DryYeast',))
    
    # Mix in Water
    M302 = units.Mixer('M302')
    P301 = units.Pump('P301')
    
    # Heat up before beer column
    # Exchange heat with stillage
    H302 = units.HXprocess('H302', outs=('', 'stillage'),
                          phase0='l', phase1='l', U=1.28)
    
    # Beer column
    x_bot = 3.910570816782338e-06
    y_top = 0.34508430224337167
    D302 = units.BinaryDistillation('D302', P=101325,
                                y_top=y_top, x_bot=x_bot, k=1.25,
                                LHK=('Ethanol', 'Water'))
    D302.tray_material = 'Stainless steel 304'
    D302.vessel_material = 'Stainless steel 304'
    D302.boiler.U = 1.85
    P302 = units.Pump('P302')
    
    # Mix ethanol Recycle (Set-up)
    M303 = units.Mixer('M303')
    x_bot = 3.910570816782338e-06
    y_top = 0.7906528373264998
    D303 = units.BinaryDistillation('D303', P=101325,
                                y_top=y_top, x_bot=x_bot, k=1.25,
                                LHK=('Ethanol', 'Water'),
                                tray_material='Stainless steel 304',
                                vessel_material='Stainless steel 304',
                                is_divided=True)
    D303.boiler.U = 1.85
    P303 = units.Pump('P303')
    
    # Superheat vapor for mol sieve
    H303 = units.HXutility('H303', T=115+273.15, V=1)
    
    # Molecular sieve
    U301 = units.MolecularSieve('U301',
                                split=(2165.14/13356.04, 1280.06/1383.85),
                                order=('Ethanol', 'Water'))
    
    sugar_solution-R301-1-T301-0-C301
    (C301-0, D301-1)-M302-P301
    (P301-0, P302-0)-H302-0-D302-1-P302
    (D302-0, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303