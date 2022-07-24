# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

__all__ = (
    'create_facilities', 
    'create_all_facilities',
    'create_coheat_and_power_system',
)

def create_all_facilities(
    feedstock,
    CIP_over_feedstock=None,
    plant_air_over_feedstock=None,
    fire_water_over_feedstock=None,
    recycle_process_water_streams=None,
    treated_water_streams=None,
    CT=True,
    CWP=True,
    CIP=True,
    FWT=True,
    ADP=True,
    WWT=True,
    CHP=True,
    HXN=True,
    PWC=True,
    CT_kwargs=None,
    CWP_kwargs=None,
    CIP_kwargs=None,
    FWT_kwargs=None,
    ADP_kwargs=None,
    WWT_kwargs=None,
    CHP_kwargs=None,
    HXN_kwargs=None,
    PWC_kwargs=None,
    area=None,
    blowdown_recycle=False,
):
    return create_facilities(
        feedstock, CIP_over_feedstock, plant_air_over_feedstock, fire_water_over_feedstock,
        recycle_process_water_streams, treated_water_streams,
        CT=CT,
        CWP=CWP,
        CIP=CIP,
        FWT=FWT,
        ADP=ADP,
        WWT=WWT,
        CHP=CHP,
        HXN=HXN,
        PWC=PWC,
        CT_kwargs=CT_kwargs,
        CWP_kwargs=CWP_kwargs,
        CIP_kwargs=CIP_kwargs,
        FWT_kwargs=FWT_kwargs,
        ADP_kwargs=ADP_kwargs,
        WWT_kwargs=WWT_kwargs,
        CHP_kwargs=CHP_kwargs,
        HXN_kwargs=HXN_kwargs,
        PWC_kwargs=PWC_kwargs,
        area=area,
        blowdown_recycle=blowdown_recycle,
    )

def create_facilities(
        feedstock=None,
        CIP_over_feedstock=None,
        plant_air_over_feedstock=None,
        fire_water_over_feedstock=None,
        recycle_process_water_streams=None,
        treated_water_streams=None,
        CT=None,
        CWP=None,
        CIP=None,
        FWT=None,
        ADP=None,
        WWT=None,
        CHP=None,
        HXN=None,
        PWC=None,
        CT_kwargs=None,
        CWP_kwargs=None,
        CIP_kwargs=None,
        FWT_kwargs=None,
        ADP_kwargs=None,
        WWT_kwargs=None,
        CHP_kwargs=None,
        HXN_kwargs=None,
        PWC_kwargs=None,
        area=None,
        blowdown_recycle=False,
    ):
    kwargs = (CHP_kwargs, WWT_kwargs, CT_kwargs, CWP_kwargs, CIP_kwargs,
              FWT_kwargs, ADP_kwargs, HXN_kwargs, PWC_kwargs)
    (CHP_kwargs, WWT_kwargs, CT_kwargs, CWP_kwargs, CIP_kwargs,
     FWT_kwargs, ADP_kwargs, HXN_kwargs, PWC_kwargs) = kwargs = [
         {} if i is None else i for i in kwargs
    ]
    if area:
        for i in kwargs[:2]: 
            if 'area' not in i: i['area'] = area
        for i in kwargs[2:]: 
            if 'ID' not in i: i['ID'] = area
    if CT: bst.facilities.CoolingTower(**CT_kwargs)
    if CWP: bst.facilities.ChilledWaterPackage(**CWP_kwargs)
    if WWT: WWT = bst.create_wastewater_treatment_system(mockup=True, autopopulate=True, **WWT_kwargs)
    if HXN: bst.HeatExchangerNetwork(**HXN_kwargs)
    if recycle_process_water_streams is None: 
        recycle_process_water_streams = ()
    if feedstock is not None:
        if CIP:
            CIP = bst.Stream('CIP', Water=126, units='kg/hr')
            CIP_package = bst.CIPpackage(ins=CIP, **CIP_kwargs)
            CIP_package.CIP_over_feedstock = 0.00121 if CIP_over_feedstock is None else CIP_over_feedstock
            @CIP_package.add_specification(run=True)
            def adjust_CIP(): CIP.imass['Water'] = feedstock.F_mass * CIP_package.CIP_over_feedstock
        if ADP:
            plant_air = bst.Stream('plant_air', N2=83333, units='kg/hr')
            ADP = bst.AirDistributionPackage(ins=plant_air, **ADP_kwargs)
            ADP.plant_air_over_feedstock = 0.8 if plant_air_over_feedstock is None else plant_air_over_feedstock
            @ADP.add_specification(run=True)
            def adjust_plant_air(): plant_air.imass['N2'] = feedstock.F_mass * ADP.plant_air_over_feedstock
        if FWT:
            fire_water = bst.Stream('fire_water', Water=8343, units='kg/hr')
            FT = bst.FireWaterTank(ins=fire_water, **FWT_kwargs)
            FT.fire_water_over_feedstock = 0.08 if fire_water_over_feedstock is None else fire_water_over_feedstock
            @FT.add_specification(run=True)
            def adjust_fire_water(): fire_water.imass['Water'] = feedstock.F_mass * FT.fire_water_over_feedstock
    if CHP:
        create_coheat_and_power_system(mockup=True, autopopulate=True, **CHP_kwargs)
    if blowdown_recycle:
        units = bst.main_flowsheet.unit.get_context_level(0)
        blowdown_to_wastewater = bst.Stream('blowdown_to_wastewater')
        bst.BlowdownMixer(
            area or '',
            [i.blowdown_water for i in units if hasattr(i, 'blowdown_water')],
            blowdown_to_wastewater
        )
    if PWC:
        process_water_mixer = bst.Mixer(area or '', ins=recycle_process_water_streams or '')
        process_water = process_water_mixer.outs[0]
        if treated_water_streams is None:
            units = bst.main_flowsheet.unit.get_context_level(0)
            treated_water_streams = [i.treated_water for i in units if hasattr(i, 'treated_water')]
        treated_water_mixer = bst.Mixer(area or '', ins=treated_water_streams)
        treated_water = treated_water_mixer.outs[0]
        bst.ProcessWaterCenter(area or '', ins=[treated_water, '', process_water], **PWC_kwargs)

@bst.SystemFactory(
    ID='CHP_sys',
    ins=['makeup_water', 'natural_gas', 'FGD_lime', 'boiler_chemicals'],
    outs=['emissions', 'blowdown', 'ash'],
    fixed_ins_size=False,
)
def create_coheat_and_power_system(
        ins, outs, combustible_slurries=None, combustible_gases=None, 
        autopopulate=None, **kwargs
    ):
    makeup_water, natural_gas, lime, boiler_chemicals, = ins
    
    slurry_mixer = bst.Mixer('slurry_mixer', ins=combustible_slurries or [])
    gas_mixer = bst.Mixer('gas_mixer', ins=combustible_gases or [])
    BT = bst.BoilerTurbogenerator(
        ins=[slurry_mixer-0, gas_mixer-0, 
             makeup_water, natural_gas, lime, boiler_chemicals],
        outs=outs,
        **kwargs,
    )
    BT.autopopulate = False if autopopulate is None else autopopulate
    
    @BT.add_specification(run=True)
    def autopopulate_combustibles():
        if BT.autopopulate and not slurry_mixer.ins and not gas_mixer.ins:
            streams = bst.FreeProductStreams(BT.system.streams)
            slurry_mixer.ins.extend(streams.combustible_slurries)
            gas_mixer.ins.extend(streams.combustible_gases)
            slurry_mixer.run_until(BT)
            gas_mixer.run_until(BT)
            BT._system.update_configuration()
