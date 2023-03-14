#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import math
import biosteam as bst
from biosteam.units.design_tools.mechanical import (
    brake_efficiency as brake_eff,
    motor_efficiency as motor_eff
    )
from . import select_pipe, format_str

__all__ = ('WWTpump',)

_hp_to_kW = 0.7457 # auom('hp').conversion_factor('kW')
_lb_to_kg = 0.4536 # auom('lb').conversion_factor('kg')
_ft_to_m = 0.3048 # auom('ft').conversion_factor('m')
_ft3_to_gal = 7.4805 # auom('ft3').conversion_factor('gallon')
_m3_to_gal = 264.1721 # auom('m3').conversion_factor('gallon')


# %%

class WWTpump(bst.Unit):
    '''
    Generic class for pumps used in wastewater treatment.

    Parameters
    ----------
    pump_type : str
        The type of the pump that determines the design algorithms to use.
        The following combination is valid:

            - "permeate_cross-flow"
            - "retentate_CSTR"
            - "retentate_AF"
            - "recirculation_CSTR"
            - "recirculation_AF"
            - "lift"
            - "sludge"
            - "chemical"

    Q_mgd : float
        Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
    add_inputs : dict
        Additional inputs that will be passed to the corresponding design algorithm.
        Check the document for the design algorithm for the specific input requirements.

    References
    ----------
    .. [1] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
        Valorization of Dilute Organic Carbon Waste Streams.
        Energy Environ. Sci. 2016, 9 (3), 1102–1112.
        https://doi.org/10.1039/C5EE03715H.
    '''
    _N_ins = 1
    _N_outs = 1

    v = 3 # fluid velocity, [ft/s]
    C = 110 # Hazen- Williams coefficient for stainless steel (SS)

    _default_equipment_lifetime = {'Pump': 15}

    _valid_pump_types = (
        'permeate_cross-flow',
        'retentate_CSTR',
        'retentate_AF',
        'recirculation_CSTR',
        'recirculation_AF',
        'lift',
        'sludge',
        'chemical'
        )


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 pump_type, Q_mgd=None, add_inputs):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.pump_type = pump_type
        self.Q_mgd = Q_mgd
        self.add_inputs = add_inputs


    def _run(self):
        self.outs[0].copy_like(self.ins[0])


    def _design(self):
        pump_type = format_str(self.pump_type)
        design_func = getattr(self, f'design_{pump_type}')

        D = self.design_results
        pipe, pumps, hdpe = design_func()
        D['Pipe stainless steel [kg]'] = pipe
        #!!! need to consider pump's lifetime in LCA
        D['Pump stainless steel [kg]'] = pumps
        D['Chemical storage HDPE [m3]'] = hdpe


    def _cost(self):
        self.power_utility.rate = self.BHP/self.motor_efficiency * _hp_to_kW


    # Used by other classes
    @staticmethod
    def _batch_adding_pump(obj, IDs, ins_dct, type_dct, inputs_dct):
        for i in IDs:
            if not hasattr(obj, f'{i}_pump'):
                if getattr(obj, 'system', None):                    
                    if not bst.main_flowsheet is obj.system.flowsheet:
                        bst.main_flowsheet.set_flowsheet(obj.system.flowsheet)
                # Add '.' in ID for auxiliary units
                pump = WWTpump(
                    ID=f'.{obj.ID}_{i}',
                    ins=ins_dct[i],
                    pump_type=type_dct[i],
                    add_inputs=inputs_dct[i])
                setattr(obj, f'{i}_pump', pump)


    # Generic algorithms that will be called by all design functions
    def _design_generic(self, Q_mgd, N_pump, L_s, L_d, H_ts, H_p):
        self.Q_mgd, self._H_ts, self._H_p = Q_mgd, H_ts, H_p
        v, C, Q_cfs = self.v, self.C, self.Q_cfs # [ft/s], -, [ft3/s]

        ### Suction side ###
        # Suction pipe (permeate header) dimensions
        OD_s, t_s, ID_s = select_pipe(Q_cfs/N_pump, v) # [in]

        # Suction friction head, [ft]
        self._H_sf = 3.02 * L_s * (v**1.85) * (C**(-1.85)) * ((ID_s/12)**(-1.17))

        ### Discharge side ###
        # Discharge pipe (permeate collector) dimensions
        OD_d, t_d, ID_d = select_pipe(Q_cfs, v)

        # Discharge friction head, [ft]
        self._H_df = 3.02 * L_d * (v**1.85) * (C**(-1.85)) * ((ID_d/12)**(-1.17))

        ### Material usage ###
        # Pipe SS, assume stainless steel, density = 0.29 lbs/in3
        # SS volume for suction, [in3]
        self._N_pump = N_pump
        V_s = N_pump * math.pi/4*((OD_s)**2-(ID_s)**2) * (L_s*12)
        # SS volume for discharge, [in3]
        V_d =math.pi/4*((OD_d)**2-(ID_d)**2) * (L_d*12)
        # Total SS mass, [kg]
        M_SS_pipe = 0.29 * (V_s+V_d) * _lb_to_kg

        # Pump SS (for pumps within 300-1000 gpm)
        # https://www.godwinpumps.com/images/uploads/ProductCatalog_Nov_2011_spread2.pdf
        # assume 50% of the product weight is SS
        M_SS_pump = N_pump * (725*0.5)

        return M_SS_pipe, M_SS_pump


    def design_permeate_cross_flow(self, Q_mgd=None, cas_per_tank=None, D_tank=None,
                                   TMP=None, include_aerobic_filter=False):
        '''
        Design pump for the permeate stream of cross-flow membrane configuration.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        cas_per_tank : int
            Number of membrane cassettes per tank.
        D_tank: float
            Depth of the membrane tank, [ft].
        TMP : float
            Transmembrane pressure, [psi].
        include_aerobic_filter : bool
            Whether aerobic filter is included in the reactor design.
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        cas_per_tank = cas_per_tank or add_inputs[0]
        D_tank = D_tank or add_inputs[1]
        TMP = TMP or add_inputs[2]
        include_aerobic_filter = include_aerobic_filter or add_inputs[3]

        H_ts_PERM = D_tank if include_aerobic_filter else 0

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=cas_per_tank,
            L_s=20, # based on a 30-module unit with a total length of 6 m, [ft]
            L_d=10*cas_per_tank, # based on a 30-module unit with a total width of 1.6 m and extra space, [ft]
            H_ts=H_ts_PERM, #  H_ds_PERM (D_tank) - H_ss_PERM (0 or D_tank)
            H_p=TMP*2.31 # TMP in water head, [ft], comment below on 2.31
            )

        # # factor = 2.31 calculated by
        # factor = auom('psi').conversion_factor('Pa') # Pa is kg/m/s2, now in [Pa]
        # factor /= 9.81 # divided by the standard gravity in m/s2, now in [kg/m2]
        # factor /= 1e3 # divided by water's density in kg/m3, now in [m]
        # factor *= auom('m').conversion_factor('ft') # m to ft

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_retentate_CSTR(self, Q_mgd=None, cas_per_tank=None):
        '''
        Design pump for the retentate stream of CSTR reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        cas_per_tank : int
            Number of membrane cassettes per tank.
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        cas_per_tank = cas_per_tank or self.add_inputs[0]

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=cas_per_tank,
            L_s=100, # pipe length per module
            L_d=30, # pipe length per module (same as the discharge side of lift pump)
            H_ts=0., # H_ds_IR (D_tank) - H_ss_IR (D_tank)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_retentate_AF(self, Q_mgd=None, N_filter=None, D=None):
        '''
        Design pump for the retentate stream of AF reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_filter : float
            Number of filter tanks.
        D : float
            Depth of the filter tank, [ft].
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_filter = N_filter or add_inputs[0]
        D = D or add_inputs[1]

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=N_filter,
            L_s=100, # assumed pipe length per filter, [ft]
            L_d=30, # same as discharge side of lift pumping, [ft]
            H_ts=0., # H_ds_IR (D) - H_ss_IR (D)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_recirculation_CSTR(self, Q_mgd=None, L_CSTR=None):
        '''
        Design pump for the recirculation stream of CSTR reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        L_CSTR : float
            Length of the CSTR tank, [ft].
        '''
        Q_mgd = Q_mgd or self.Q_mgd
        L_CSTR = L_CSTR or self.add_inputs[0]

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=1,
            L_s=0., # ignore suction side
            L_d=L_CSTR, # pipe length per train
            H_ts=5., # H_ds_IR (5) - H_ss_IR (0)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_recirculation_AF(self, Q_mgd=None, N_filter=None, d=None,
                                D=None):
        '''
        Design pump for the recirculation stream of AF reactors.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_filter : float
            Number of filter tanks.
        d : float
            diameter of the filter tank, [ft].
        D : float
            Depth of the filter tank, [ft].
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_filter = N_filter or add_inputs[0]
        d = d or add_inputs[1]
        D = D or add_inputs[2]

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=N_filter,
            L_s=d+D, # pipe length per filter, [ft]
            L_d=30, # same as discharge side of lift pumping, [ft]
            H_ts=0., # H_ds_IR (D) - H_ss_IR (D)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_lift(self, Q_mgd=None, N_filter=None, D=None):
        '''
        Design pump for the filter tank to lift streams.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        N_filter : float
            Number of filter tanks.
        D : float
            Depth of the filter tank, [ft].
        '''
        add_inputs = self.add_inputs
        Q_mgd = Q_mgd or self.Q_mgd
        N_filter = N_filter or add_inputs[0]
        D = D or add_inputs[1]

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=N_filter,
            L_s=150, # length of suction pipe per filter, [ft]
            L_d=30, # pipe length per filter, [ft]
            H_ts=D, # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_sludge(self, Q_mgd=None):
        '''
        Design pump for handling waste sludge.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        '''
        Q_mgd = Q_mgd or self.Q_mgd

        M_SS_IR_pipe, M_SS_IR_pump = self._design_generic(
            Q_mgd=Q_mgd,
            N_pump=1,
            L_s=50, # length of suction pipe, [ft]
            L_d=50, # length of discharge pipe, [ft]
            H_ts=0., # H_ds_LIFT (D) - H_ss_LIFT (0)
            H_p=0. # no pressure
            )

        return M_SS_IR_pipe, M_SS_IR_pump, 0


    def design_chemical(self, Q_mgd=None):
        '''
        Design pump for membrane cleaning chemicals (NaOCl and citric acid),
        storage containers are included, and are assumed to be cubic in shape
        and made of HDPE.

        Parameters defined through the `add_inputs` argument upon initialization of
        this unit (Q_mgd listed separately) will be used if not provided
        when calling this function.

        Parameters
        ----------
        Q_mgd : float
            Volumetric flow rate in million gallon per day, [mgd].
        '''
        if not Q_mgd:
            V_CHEM = self.ins[0].F_vol * 24 * 7 * 2 # for two weeks of storage, [m3]
            Q_CHEM_mgd = self.Q_mgd
        else:
            V_CHEM = (Q_mgd*1e6/_m3_to_gal) * 7 * 2
            Q_CHEM_mgd = Q_mgd

        # HDPE volume, [m3], 0.003 [m] is the thickness of the container
        V_HDPE = 0.003 * (V_CHEM**(1/3))**2*6
        # # Mass of HDPE, [m3], 950 is the density of the HDPE in [kg/m3]
        # M_HDPE = 950 * V_HDPE

        H_ss_CHEM = V_CHEM**(1/3) / _ft_to_m
        # 9'-7" is the water level in membrane trains
        # 18" is the distance from C/L of the pump to the ground
        H_ds_CHEM = 9 + 7/12 - 18/12
        H_ts_CHEM = H_ds_CHEM - H_ss_CHEM

        M_SS_CHEM_pipe, M_SS_CHEM_pump = self._design_generic(
            Q_mgd=Q_CHEM_mgd,
            N_pump=1,
            L_s=0., # no suction pipe
            L_d=30.,
            H_ts=H_ts_CHEM,
            H_p=0. # no pressure
            )

        return M_SS_CHEM_pipe, M_SS_CHEM_pump, V_HDPE


    @property
    def pump_type(self):
        '''
        [str] The type of the pump that determines the design algorithms to use.
        Use `valid_pump_type` to see acceptable pump types.
        '''
        return self._pump_type
    @pump_type.setter
    def pump_type(self, i):
        i_lower = i.lower()
        i_lower = i_lower.replace('cstr', 'CSTR')
        i_lower = i_lower.replace('af', 'AF')
        if i_lower not in self.valid_pump_types:
            raise ValueError(f'The given `pump_type` "{i}" is not valid, '
                             'check `valid_pump_types` for acceptable pump types.')
        self._pump_type = i_lower

    @property
    def valid_pump_types(self):
        '''[tuple] Acceptable pump types.'''
        return self._valid_pump_types

    @property
    def N_pump(self):
        '''[int] Number of pumps.'''
        return self._N_pump

    @property
    def H_sf(self):
        '''[float] Suction friction head, [ft].'''
        return self._H_sf

    @property
    def H_df(self):
        '''[float] Discharge friction head, [ft].'''
        return self._H_df

    @property
    def H_ts(self):
        '''[float] Total static head, [ft].'''
        return self._H_ts

    @property
    def H_p(self):
        '''[float] Pressure head, [ft].'''
        return self._H_p

    @property
    def TDH(self):
        '''[float] Total dynamic head, [ft].'''
        return self.H_ts+self.H_sf+self.H_df+self.H_p

    @property
    def BHP(self):
        '''[float] Brake horsepower, [hp].'''
        return (self.TDH*self.Q_gpm)/3960/self.brake_efficiency

    @property
    def Q_mgd(self):
        '''
        [float] Volumetric flow rate in million gallon per day, [mgd].
        Will use total volumetric flow through the unit if not provided.
        '''
        if self._Q_mgd:
            return self._Q_mgd
        return self.F_vol_in*_m3_to_gal*24/1e6
    @Q_mgd.setter
    def Q_mgd(self, i):
        self._Q_mgd = i

    @property
    def Q_gpm(self):
        '''[float] Volumetric flow rate in gallon per minute, [gpm].'''
        return self.Q_mgd*1e6/24/60

    @property
    def Q_cmd(self):
        '''
        [float] Volumetric flow rate in cubic meter per day, [cmd].
        '''
        return self.Q_mgd *1e6/_m3_to_gal # [m3/day]

    @property
    def Q_cfs(self):
        '''[float] Volumetric flow rate in cubic feet per second, [cfs].'''
        return self.Q_mgd*1e6/24/60/60/_ft3_to_gal

    @property
    def brake_efficiency(self):
        '''[float] Brake efficiency.'''
        return brake_eff(self.Q_gpm)

    @property
    def motor_efficiency(self):
        '''[float] Motor efficiency.'''
        return motor_eff(self.BHP)