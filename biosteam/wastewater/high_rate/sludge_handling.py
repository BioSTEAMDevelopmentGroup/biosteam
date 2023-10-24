# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import math, flexsolve as flx
import biosteam as bst
from warnings import warn

__all__ = ('SludgeHandling', 'BeltThickener', 'SludgeCentrifuge')


# %%

class SludgeHandling(bst.Unit):
    """
    A generic class for handling of wastewater treatment sludge [6]_. Two pumps 
    (one for the supernatant and one for sludge) are included.

    Separation split is determined by the moisture (i.e., water) content
    of the sludge, soluble chemicals will have the same split as water,
    insolubles chemicals will all go to the retentate.

    Parameters
    ----------
    ins : 
        Influent
    outs : 
        * [0] water-rich supernatant
        * [1] solid-rich sludge
    sludge_moisture : float
        Moisture content of the sludge, [wt. % water].
    solubles : tuple
        IDs of the soluble chemicals.
        Note that all chemicals that are not included in this tuple and not
        locked as gas phase (i.e., `chemical.locked_state` != "g") will be
        treated as solids in simulation.

    """

    SKIPPED = False
    _graphics = bst.Splitter._graphics
    _ins_size_is_fixed = False
    _N_outs = 2
    auxiliary_unit_names = ('effluent_pump', 'sludge_pump')


    def _init(self, sludge_moisture=0.96, solubles=()):
        self.sludge_moisture = sludge_moisture
        self.solubles = tuple(solubles)
        self.solids = tuple(i.ID for i in self.chemicals
                            if (i.ID not in solubles) and (i.locked_state!='g'))
        ID = self.ID
        self._mixed = bst.Stream(f'{ID}_mixed')
        # Add '.' in ID for auxiliary units
        self.effluent_pump = bst.Pump(f'.{ID}_eff_pump', ins=self.outs[0].proxy(f'{ID}_eff'))
        self.sludge_pump = bst.Pump(f'.{ID}_sludge_pump', ins=self.outs[1].proxy(f'{ID}_sludge'))


    @staticmethod
    def _mc_at_split(split, solubles, mixed, eff, sludge, target_mc):
        eff.imass[solubles] = mixed.imass[solubles] * split
        sludge.imass[solubles] = mixed.imass[solubles] - eff.imass[solubles]
        mc = sludge.imass['Water'] / sludge.F_mass
        return mc-target_mc


    def _run(self):
        eff, sludge = self.outs
        solubles, solids = self.solubles, self.solids

        mixed = self._mixed
        mixed.mix_from(self.ins)
        eff.T = sludge.T = mixed.T

        if mixed.F_mass == 0:
            eff.empty()
            sludge.empty()
            self.SKIPPED = True
            return

        mc = mixed.imass['Water']/mixed.F_mass
        target_mc = self.sludge_moisture
        if mc < target_mc:
            if self.SKIPPED == False:
                warn(f'The moisture of influent ({mc:.2%}) is smaller than the desired {target_mc:.0%}, '
                     'simulation skipped.')
                self.SKIPPED = True
            sludge.copy_like(mixed)
            eff.empty()
            return

        sludge.copy_flow(mixed, solids, remove=True) # all solids go to sludge
        eff.copy_flow(mixed, solubles)

        flx.IQ_interpolation(
                f=self._mc_at_split, x0=1e-3, x1=1.-1e-3, ytol=1e-3, maxiter=100,
                args=(solubles, mixed, eff, sludge, target_mc),
                checkbounds=False)
        self.SKIPPED = False


    def _cost(self):
        if self.SKIPPED == False:
            pumps = (self.effluent_pump, self.sludge_pump)
            for i in range(2): pumps[i].simulate()
        else:
            self.baseline_purchase_costs.clear()
            self.power_utility.rate = 0


class BeltThickener(SludgeHandling):
    """
    Gravity belt thickener (GBT) designed based on the manufacture specification
    data sheet. [9]_

    Key parameters include:
    
    * Capacity: 80-100 m3/h.
    * Influent solids concentration: 0.2-1%.
    * Sludge cake moisture content: 90-96%.
    * Motor power: 3 (driving motor) and 1.1 (agitator motor) kW.
    * Belt width: 2.5 m.
    * Weight: 2350 kg.
    * Quote price: $3680 ea for three or more sets.

    The bare module (installation) factor is from Table 25 in Humbird et al. [2]_
    (solids handling equipment).

    Parameters
    ----------
    ins : 
        Influent
    outs : 
        * [0] water-rich supernatant
        * [1] solid-rich sludge
    sludge_moisture : float
        Moisture content of the thickened sludge, [wt% water].
    solubles : tuple
        IDs of the soluble chemicals.
        Note that all chemicals that are not included in this tuple and not
        locked as gas phase (i.e., `chemical.locked_state!='g'`) will be
        treated as solids in simulation.
    max_capacity : float
        Maximum hydraulic loading per belt thickener, [m3/h].
    power_demand : float
        Total power demand of each belt thickener, [kW].

    """

    def _init(self, sludge_moisture=0.96, solubles=(),
                 max_capacity=100, power_demand=4.1):
        SludgeHandling._init(self, sludge_moisture=sludge_moisture,
                             solubles=solubles)
        self.max_capacity = max_capacity
        self.power_demand = power_demand


    def _design(self):
        self._N_thickener = N = math.ceil(self._mixed.F_vol/self.max_capacity)
        self.design_results['Number of thickeners'] = N
        self.F_BM['Thickeners'] = 1.7 # ref [2]
        self.baseline_purchase_costs['Thickeners'] = 4000 * N


    def _cost(self):
        super()._cost()
        self.power_utility.rate += self.power_demand * self.N_thickener


    @property
    def N_thickener(self):
        """[int] Number of required belt thickeners."""
        return self._N_thickener


class SludgeCentrifuge(SludgeHandling, bst.SolidsCentrifuge):
    """
    Solid centrifuge for sludge dewatering.

    The simulation (:meth:`~biosteam.Unit._run` method) and costing (:meth:`~biosteam.Unit._cost` method) are based on 
    :class:`~biosteam.SludgeHandling` and the sizing (:meth:`~biosteam.Unit._design` method) is based on 
    :class:`~biosteam.units.solids_separation.SolidsCentrifuge`.

    Parameters
    ----------
    ins : 
        Influent
    outs : 
        * [0] water-rich supernatant
        * [1] solid-rich sludge
    sludge_moisture : float
        Moisture content of the thickened sludge, [wt% water].
    solubles : tuple
        IDs of the soluble chemicals.
        Note that all chemicals that are not included in this tuple and not
        locked as gas phase (i.e., `chemical.locked_state!='g'`) will be
        treated as solids in simulation.
    
    """

    def _init(self, sludge_moisture=0.8, solubles=(),
              centrifuge_type='scroll_solid_bowl'):
        SludgeHandling._init(self, sludge_moisture=sludge_moisture,
                             solubles=solubles)
        self.centrifuge_type = centrifuge_type

    _run = SludgeHandling._run

    _design = bst.SolidsCentrifuge._design

    _cost = SludgeHandling._cost