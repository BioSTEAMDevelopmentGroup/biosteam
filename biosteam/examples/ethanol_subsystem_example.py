# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biosteam import units
from thermosteam import functional as fn
from biosteam import main_flowsheet
import biosteam as bst
import thermosteam as tmo

__all__ = ('ethanol_subsystem_example',
)

def ethanol_subsystem_example():
    """
    Test BioSTEAM by creating a conventional sugarcane fermentation and ethanol
    purification process.
    
    Examples
    --------
    >>> ethanol_sys = ethanol_subsystem_example()
    >>> # The sugarcane_example_subsystem flowsheet may help for accessing units
    >>> from biosteam import main_flowsheet as F
    >>> fs = F.flowsheet['ethanol_subsystem_example']
    >>> fs.unit # Check unit operation registry
    Register:
     <Fermentation: R301>
     <StorageTank: T301>
     <VentScrubber: D301>
     <SolidsCentrifuge: C301>
     <Mixer: M302>
     <Pump: P301>
     <HXprocess: H302>
     <BinaryDistillation: D302>
     <Pump: P302>
     <Mixer: M303>
     <BinaryDistillation: D303>
     <Pump: P303>
     <HXutility: H303>
     <MolecularSieve: U301>
    >>> R301 = fs.unit.R301 # Get unit operation

    """
    original_flowsheet = main_flowsheet.get_flowsheet()
    main_flowsheet.set_flowsheet('ethanol_subsystem_example')
    
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
    fermentation_feed = bst.Stream('fermentation_feed', phase='l',
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
    
    fermentation_feed-R301-1-T301-0-C301
    (C301-0, D301-1)-M302-P301
    (P301-0, P302-0)-H302-0-D302-1-P302
    (D302-0, U301-0)-M303-0-D303-0-H303-U301
    D303-1-P303
    
    ethanol_subsystem_example = main_flowsheet.create_system('ethanol_subsystem_example')
    ethanol_subsystem_example.simulate()
    main_flowsheet.set_flowsheet(original_flowsheet) 
    return ethanol_subsystem_example
    