# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-2024, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
sympy = None
import biosteam as bst
from . import (
    get_BD_dct,
    compute_stream_COD,
    get_digestion_rxns,
    IC_purchase_cost_algorithms
)


__all__ = ('InternalCirculationRx',)


# %%

class InternalCirculationRx(bst.MixTank):
    """
    Internal circulation (IC) reactor for anaerobic digestion (AD),
    including a high-rate bottom reactor for rapid organic removal and
    a low-rate top reactor for polishing.
    Both reactors are similar to upflow anaerobic blanket reactor (UASB).

    Design of the reactor follows steps described in [4]_
    (assuming steady state and pseudo-zeroth-order kinetics),
    where two methods are used based on Irizar et al. [7]_ and
    Tchobanoglous et al. [8]_.

    Parameters
    ----------
    ins : 
        Influent.
    outs : 
        * [0] biogas
        * [1] effluent
        * [2] waste sludge
    method : str
        * "separate" to design the bottom and top reactors separately as in [7]_.
        Design parameters for this methid include OLRall, biodegradability, Y, 
        q_Qw, mu_max, b, Fxt, and Fxb.
        
        * "lumped" to design the entire IC reactor as a black box following [8]_.
        Design parameters for this method include OLRall, biodegradability, Y,
        q_Qw, and q_Xw.
    OLRall : float
        Overall organic loading rate, [kg COD/m3/hr].
    Y_biogas : float
        Biogas yield, [kg biogas/kg consumed COD].
    Y_biomass : float
        Biomass yield, [kg biomass/kg consumed COD].
    biodegradability : float or dict
        Biodegradability of chemicals,
        when shown as a float, all biodegradable chemicals are assumed to have
        the same degradability.
    q_Qw : float
        Ratio between the bottom reactor waste flow and the influent.
    q_Xw : float
        Ratio between the biomass concentration in the reactor and the waste flow.
    mu_max : float
        Maximum specific growth rate, [/hr].
    b : float
        Specific endogenous decay coefficient, [/hr].
    V_wf : float
        Fraction of working volume over total volume.
    vessel_type : str
        Can be "IC" to use the reactor size constraints according to [4]_,
        or "Conventional" based on :class:`biosteam.MixTank`
        (much smaller tank size, not recommended).
    vessel_material : str
        Vessel material.
    kW_per_m3 : float
        Electricity requirement per unit volume, [kW/m^3].
        Default to 0 as IC reactors realizes mixing through internal circulation
        caused by the rising force of the generated biogas.
    T : float
        Temperature of the reactor.
        Will not control temperature if provided as None.
    kwargs : dict
        Other keyword arguments (e.g., Fxb, Fxt).

    """
    _N_ins = 1
    _N_outs = 3 # biogas, effluent, waste sludge

    # Assumptions
    _q_Qw = 0.01
    _q_Xw = 1.5
    _mu_max = 0.01
    _b = 0.00083
    _Fxb = 0.0032
    _Fxt = 0.0281
    _Y = 0.05

    # Related to cost algorithm
    _default_vessel_type = 'IC'
    _default_vessel_material = 'Stainless steel'
    purchase_cost_algorithms = IC_purchase_cost_algorithms

    # Other equipment
    auxiliary_unit_names = ('heat_exchanger', 'effluent_pump', 'sludge_pump')


    def _init(self, method='lumped', OLRall=1.25, Y_biogas=0.86, Y_biomass=0.05, biodegradability={},
              vessel_type='IC', vessel_material='Stainless steel',
              V_wf=0.8, kW_per_m3=0., T=35+273.15, hxn_ok=False, **kwargs):
        self.method = method
        self.OLRall = OLRall
        self.Y_biogas = Y_biogas
        self.Y_biomass = Y_biomass
        self.biodegradability = \
            biodegradability if biodegradability else get_BD_dct(self.chemicals)
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_type = 'IC'
        self.vessel_material = vessel_material
        self.kW_per_m3 = kW_per_m3
        self.T = T
        # Initialize the attributes
        ID = self.ID
        self._inf = bst.Stream(f'{ID}_inf')
        hx_in = bst.Stream(f'{ID}_hx_in')
        hx_out = bst.Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = bst.HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)
        self._refresh_rxns()
        # Conversion will be adjusted in the _run function
        self._decay_rxn = self.chemicals.WWTsludge.get_combustion_reaction(conversion=0.)
        self.effluent_pump = bst.Pump(f'.{ID}_eff', ins=self.outs[1].proxy(f'{ID}_eff'))
        self.sludge_pump = bst.Pump(f'.{ID}_sludge', ins=self.outs[2].proxy(f'{ID}_sludge'))
        self.hxn_ok = hxn_ok
        for k, v in kwargs.items(): setattr(self, k, v)

    def _refresh_rxns(self, Y_biogas=None, Y_biomass=None):
        Y_biogas = Y_biogas if Y_biogas else self.Y_biogas
        Y_biomass = Y_biomass if Y_biomass else self.Y_biomass

        self._biogas_rxns = get_digestion_rxns(self.ins[0], self.biodegradability,
                                               Y_biogas, 0., 'WWTsludge')
        self._growth_rxns = get_digestion_rxns(self.ins[0], self.biodegradability,
                                               0., Y_biomass, 'WWTsludge')
        self._i_rm = self._biogas_rxns.X + self._growth_rxns.X


    @staticmethod
    def _degassing(original_stream, receiving_stream):
        gases = tuple(i.ID for i in original_stream.chemicals if i.locked_state=='g')
        receiving_stream.imass[gases] += original_stream.imass[gases]
        original_stream.imass[gases] = 0


    @staticmethod
    def compute_COD(stream):
        r"""
        Compute the chemical oxygen demand (COD) of a given stream in kg-O2/m3
        by summing the COD of each chemical in the stream using:

        .. math::
            COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol\ chemical}
        """
        return compute_stream_COD(stream)


    def _run(self):
        inf = self._inf
        inf.copy_like(self.ins[0])
        biogas, eff, waste  = self.outs
        degassing = self._degassing

        # Initialize the streams
        biogas.phase = 'g'
        biogas.empty()

        inf.split_to(waste, eff, self.q_Qw)
        biogas_rxns = self.biogas_rxns
        growth_rxns = self.growth_rxns

        growth_rxns(inf.mol)
        biogas_rxns(inf.mol)

        degassing(inf, biogas)
        Se = self.compute_COD(inf)

        Qi, Si, Xi, Qe, Y = self.Qi, self.Si, self.Xi, self.Qe, self.Y_biomass
        method = self.method.lower()
        if method == 'separate':
            run_inputs = (Qi, Si, Xi, Qe, Se, self.Vliq, Y,
                          self.mu_max, self.b, self.Fxb, self.Fxt)
            Xw, Xe = self._run_separate(run_inputs)

        else:
            Xe = (Qi*Xi+Qi*(Si-Se)*Y)/(Qe+(Qi-Qe)*self.q_Xw)
            Xw = Xe * self.q_Xw
            for rxns in (growth_rxns, biogas_rxns):
                rxns(waste.mol)
                rxns(eff.mol)

        degassing(eff, biogas)
        degassing(waste, biogas)

        eff.imass['WWTsludge'] = Xe*self.Qe
        waste.imass['WWTsludge'] = Xw*self.Qw
        diff = sum(i.imol['WWTsludge'] for i in self.outs) - inf.imol['WWTsludge']
        if diff > 0:
            decay_rxn = self.decay_rxn
            decay_rxn._X = min(1., diff / inf.imol['WWTsludge'])

            for i in (eff, waste):
                decay_rxn.force_reaction(i.mol)
                i.imol['O2'] = max(0, i.imol['O2'])
                degassing(i, biogas)

        if self.T: biogas.T = eff.T = waste.T = self.T


    def _run_separate(self, run_inputs):
        global sympy 
        if sympy is None: import sympy
        Qi, Si, Xi, Qe, Se, Vliq, Y, mu_max, b, Fxb, Fxt = run_inputs
        Qw = Qi - Qe
        Xb, Xe, Sb, Vb = sympy.symbols('Xb, Xe, Sb, Vb', real=True)

        # Mass balances based on biomass/substrate changes in the bottom/top rx,
        # (0 at steady state)
        biomass_b = Qi*Xi - (Qe*Xb*Fxb+Qw*Xb) + Xb*Vb*(mu_max-b)
        biomass_t = Qe*(Fxb*Xb-Fxt*Xe) + Xe*(Vliq-Vb)*(mu_max-b)
        substrate_b = Qi*(Si-Sb) - mu_max*(Xb*Vb/Y)
        substrate_t = Qe*(Sb-Se) - mu_max*((Vliq-Vb)*Xe/Y)

        parameters = (Qi, Qe, Si, Se, Vliq)
        results = sympy.solve(
            (sympy.Eq(biomass_b, 0),
             sympy.Eq(biomass_t, 0),
             sympy.Eq(substrate_b, 0),
             sympy.Eq(substrate_t, 0)), (Xb, Xe, Sb, Vb))

        Xb, Xe, Sb, Vb = self._filter_results('separate', parameters, results)

        Vt = Vliq - Vb # volume of the top rx, m3
        self._Vb, self._Vt = Vb, Vt
        return Xb, Xe


    @staticmethod
    def _filter_results(method, parameters, results):
        """Check if the solution satisfies the design constraints."""
        Qi, Qe, Si, Se, Vliq = parameters
        solutions = []
        for result in results:
            Xb, Xe, Sb, Vb = result
            Vt = Vliq - Vb
            OLRt = Qe*Sb / Vt
            OLRb = Qi*Si / Vb

            if (
                    0 <= OLRt <= OLRb and
                    0 <= Se <= Sb <= Si and
                    0 <= Xe <= Xb and
                    0 <= Vb <= Vliq
                ):
                solutions.append(result)

        if len(solutions) == 0 :
            raise bst.DesignError('No feasible design found for the given parameters.')

        elif len(solutions) >1: # find more than one solution
            Xbs = [i[1] for i in solutions]
            index = Xbs.index(min(Xbs)) # choose the one with lowest effluent biomass
            return solutions[index]


    _units = {
        'HRT': 'hr',
        'SRT': 'hr',
        'Single reactor liquid volume': 'm3',
        'Bottom reactor volume': 'm3',
        'Top reactor volume': 'm3',
        'Gas chamber volume': 'm3'
        }
    def _design(self):  
        D = self.design_results
        D['HRT'] = D['Residence time'] = self.HRT
        D['SRT'] = self.SRT
        D['Total volume'] = self.Vtot
        D['Total liquid volume'] = self.Vliq
        if self.method == 'separate':
            D['Bottom reactor volume'] = self.Vb
            D['Top reactor volume'] = self.Vt


    def _cost(self):       
        bst.MixTank._cost(self)
        
        hx = self.heat_exchanger
        ins0 = self.ins[0]
        hx.ins[0].copy_flow(ins0)
        hx.outs[0].copy_flow(ins0)
        hx.ins[0].T = ins0.T
        hx.outs[0].T = self.T
        hx.ins[0].P = hx.outs[0].P = ins0.P
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs, 
                                           scale=1. / self.parallel.get('self', 1.),
                                           vle=False,
                                           hxn_ok=self.hxn_ok)
        
        for p in (self.effluent_pump, self.sludge_pump): p.simulate()



    @property
    def method(self):
        """[str] Design method, can be "separate" or "lumped"."""
        return self._method
    @method.setter
    def method(self, i):
        if not i.lower() in ('separate', 'lumped'):
            raise ValueError('`method` can only be "separated", or "lumped", '
                             f'not "{i}".')
        self._method = i.lower()

    @property
    def OLRall(self):
        """[float] Overall organic loading rate, [kg COD/m3/hr]."""
        return self._OLRall
    @OLRall.setter
    def OLRall(self, i):
        if i < 0:
            raise ValueError('`OLRall` should be >=0, '
                             f'the input value {i} is outside the range.')
        self._OLRall = i

    @property
    def biodegradability(self):
        """
        [float of dict] Biodegradability of chemicals,
        when shown as a float, all biodegradable chemicals are assumed to have
        the same degradability.
        """
        return self._biodegradability
    @biodegradability.setter
    def biodegradability(self, i):
        if not isinstance(i, dict):
            if not 0<=i<=1:
                raise ValueError('`biodegradability` should be within [0, 1], '
                                 f'the input value {i} is outside the range.')
            self._biodegradability = dict.fromkeys(self.chemicals.IDs, i)
        else:
            for k, v in i.items():
                if not 0<=v<=1:
                    raise ValueError('`biodegradability` should be within [0, 1], '
                                     f'the input value for chemical "{k}" is '
                                     'outside the range.')
            self._biodegradability = dict.fromkeys(self.chemicals.IDs, i).update(i)
        self._refresh_rxns()

    @property
    def i_rm(self):
        """[:class:`np.array`] Removal of each chemical in this reactor."""
        return self._i_rm

    @property
    def q_Qw(self):
        """[float] Ratio between the bottom reactor waste flow and the influent."""
        return self._q_Qw
    @q_Qw.setter
    def q_Qw(self, i):
        if not 0<=i<=1:
            raise ValueError('`q_Qw` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._q_Qw = i

    @property
    def q_Xw(self):
        """
        [float] Ratio between the biomass concentration in the reactor and the waste flow,
        only relevant when the "lumped" method is used.
        """
        return self._q_Xw if self.method=='lumped' else None
    @q_Xw.setter
    def q_Xw(self, i):
        if not i>=1:
            raise ValueError('`q_Xw` should be >=1, '
                             f'the input value {i} is outside the range.')
        self._q_Xw = i

    @property
    def mu_max(self):
        """
        [float] Maximum specific growth rate, [/hr],
        only relevant when the "separate" method is used.
        """
        return self._mu_max if self.method=='separate' else None
    @mu_max.setter
    def mu_max(self, i):
        if i < 0:
            raise ValueError('`mu_max` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._mu_max = i

    @property
    def b(self):
        """
        [float] Specific endogenous decay coefficient, [/hr],
        only relevant when the "separate" method is used.
        """
        return self._b if self.method=='separate' else None
    @b.setter
    def b(self, i):
        if i < 0:
            raise ValueError('`b` should be >= 0, '
                             f'the input value {i} is outside the range.')
        self._b = i

    @property
    def Fxb(self):
        """
        [float] Biomass transfer ratio from the bottom reactor to the top reactor,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        """
        return self._Fxb if self.method=='separate' else None
    @Fxb.setter
    def Fxb(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxb` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxb = i

    @property
    def Fxt(self):
        """
        [float] Biomass transfer ratio from the top reactor to the effluent,
        should be within [0, 1] (ideal to no retention),
        only relevant when the "separate" method is used.
        """
        return self._Fxt if self.method=='separate' else None
    @Fxt.setter
    def Fxt(self, i):
        if not 0<=i<=1:
            raise ValueError('`Fxt` should be within [0, 1], '
                             f'the input value {i} is outside the range.')
        self._Fxt = i

    @property
    def Vb(self):
        """
        [float] Volume of the bottom reactor, [m3],
        only relevant when the "separate" method is used.
        """
        return self._Vb if self.method=='separate' else None

    @property
    def Vt(self):
        """
        [float] Volume of the top reactor, [m3],
        only relevant when the "separate" method is used.
        """
        return self._Vt if self.method=='separate' else None

    @property
    def Vliq(self):
        """
        [float] Total volume of the liquid, not considering gas headspace
        and `V_wf`, [m3].
        """
        return self.Qi*self.Si / self.OLRall

    @property
    def Vtot(self):
        """
        [float] Total volume considering `V_wf`, [m3].
        """
        return self.Vliq / self.V_wf

    @property
    def HRT(self):
        """[float] Hydraulic retention time [hr]."""
        return self.Vliq / self.Qi

    @property
    def tau(self):
        """
        [float] Reactor residence time, [hr]
        (same as the hydraulic retention time, HRT).
        """
        return self.HRT

    @property
    def SRT(self):
        """[float] Solid residence time [hr]."""
        if self.method == 'separate':
            return (self.Xb*self.Vb+self.Xe*self.Vt)/(self.q_Qw*self.Xb+self.Fxt*self.Qe*self.Xe)

        return  self.Xe*self.Vliq / (self.q_Qw*self.q_Xw+self.Qe*self.Xe)

    @property
    def biogas_rxns(self):
        """
        [:class:`tmo.ParallelReaction`] Biogas production reactions.
        """
        return self._biogas_rxns

    @property
    def growth_rxns(self):
        """
        [:class:`tmo.ParallelReaction`] Biomass (WWTsludge) growth reactions.
        """
        return self._growth_rxns

    @property
    def decay_rxn(self):
        """
        [:class:`tmo.Reaction`] Biomass endogenous decay.

        .. note::
            Conversion is adjusted in the _run function.
        """
        return self._decay_rxn

    @property
    def Qi(self):
        """[float] Influent volumetric flow rate, [m3/hr]."""
        return self.ins[0].F_vol

    @property
    def Qe(self):
        """[float] Effluent volumetric flow rate, [m3/hr]."""
        return self.outs[1].F_vol

    @property
    def Qw(self):
        """[float] Waste flow volumetric flow rate, [m3/hr]."""
        return self.outs[2].F_vol

    @property
    def Si(self):
        """
        [float] Influent substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        """
        return self.compute_COD(self.ins[0])

    @property
    def Se(self):
        """
        [float] Effluent substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        """
        return self.compute_COD(self.outs[1])

    @property
    def Sw(self):
        """
        [float] Waste flow substrate (i.e., biodegradable chemicals)
        concentration, [kg/m3].
        """
        return self.compute_COD(self.outs[2])

    @property
    def Xi(self):
        """[float] Influent biomass (i.e., `WWTsludge`) concentration, [kg/m3]."""
        return self.ins[0].imass['WWTsludge']/self.ins[0].F_vol

    @property
    def Xe(self):
        """[float] Effluent  biomass (i.e., `WWTsludge`) concentration, [kg/m3]."""
        return self.outs[1].imass['WWTsludge']/self.outs[1].F_vol

    @property
    def Xw(self):
        """[float] Waste flow biomass (i.e., `WWTsludge`) concentration, [kg/m3]."""
        return self.outs[2].imass['WWTsludge']/self.outs[2].F_vol

    @property
    def organic_rm(self):
        """[float] Overall organic (COD) removal rate."""
        return 1 - self.Qe*self.Se/(self.Qi*self.Si)