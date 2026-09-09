# -*- coding: utf-8 -*-
"""
.. contents:: :local:

.. autoclass:: biosteam.units.aerated_bioreactor.GasFedBioreactor

References
----------
.. [1] Benz, G. T. Optimize Power Consumption in Aerobic Fermenters. 
    Chem. Eng. Progress 2003, 99 (5), 100–103.

.. [2] Benz, G. T. Bioreactor Design for Chemical Engineers. Chem. Eng.\
    Progress 2011, 21–26.

.. [3] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani, R.,
    & Ng, M. K. (2017). Product and Process Design Principles. Wiley.

"""
import biosteam as bst
from .abstract_stirred_tank_reactor import AbstractStirredTankReactor
from .aerated_bioreactor import AeratedBioreactor
from math import pi
import numpy as np
from scipy.constants import g
import flexsolve as flx
from warnings import filterwarnings, catch_warnings
from scipy.optimize import fsolve
from biosteam.units.design_tools import aeration
# from matplotlib import pyplot as plt, rcParams, colormaps

__all__ = (
    'GasFedBioreactor', 'GFB',
)

class GasFedBioreactor(AbstractStirredTankReactor):
    """
    Create a gas-fed bioreactor which satisfies the substrate mass tranfer 
    requirement of the mass balance. The gas-fed bioreactor may include 
    multiple gas feeds. An optional, user-specified `titer` can be satisfied by varying
    the flow rates of the `controlled_feeds`. If a gas feed is varied, 
    `backward_reactions` may need to be specified. 

    The reactor is designed as a pressure vessel with a given aspect ratio and 
    residence time. A pump-heat exchanger recirculation loop can be used to satisfy 
    the duty, if any. By default, a turbine agitator is also included if the 
    power usage, `kW_per_m3`, is positive. A vacuum system is also 
    automatically added if the operating pressure is at a vacuum. 

    Parameters
    ----------
    gas_substrates : 
        Substrates within the gas phase.
    titer :
        Dictionary of substrate/titer pairs [g / L].
    backward_reactions :
        Backwards reactions to get the substrate mass transfer requirement.
    controlled_gas_feeds :
        Feeds that can be varied to meet mass transfer requirement.
    theta : 
        Fraction of gas substrate saturation in the broth. Defaults to 0.5.
    Q_consumption :
        Forced duty per gas substrate consummed [kJ/kmol].
    optimize_power :
        If true, the agitator power is solved to minimize the total power 
        requirement of both the compressor and agitator such that the 
        required oxygen transfer rate is met.
    design : 
        Bioreactor design configuration. Valid options include 'Stirred tank'
        and 'Bubble column'. Defaults to the former.
    method :
        Method to calculate the overall mass transfer coefficient, kLa. 
        Can be a name or a function. Valid method names are listed in 
        `biosteam.aeration.kLa_methods`.
        For stirred tanks, defaults to the 'Riet'. 
        For bubble columns, defaults to 'Dewes'.
    kLa_kwargs:
        Additional arguments to pass to the kLa method.
    cooler_pressure_drop :
        Pressure drop at the cooler [Pa]. Defaults to 20684.28 Pa, 
        a heuristic value for a gas.
    compressor_isentropic_efficiency :
        Isentropic efficiency of the compressor. Defaults to 0.85.
    tau :
        Residence time [hr].
    T : 
        Operating temperature [K].
    P : 
        Operating pressure [Pa].
    V_wf : 
        Fraction of working volume over total volume. Defaults to 0.8.
    length_to_diameter :
        Length to diameter ratio of bioreactor.
    V_max :
        Maximum volume of a reactor [m3]. Defaults to 355.
    kW_per_m3 : 
        Power usage of agitator. Defaults to 0.985 [kW / m3] converted from 
        5 hp/1000 gal as in [1]_, for liquid–liquid reaction or extraction.
    vessel_material : 
        Vessel material. Defaults to 'Stainless steel 316'.
    vessel_type : 
        Vessel type. Valid options are 'Horizontal' or 'Vertical'. Defaults to 'Vertical'
    batch :
        Whether to use batch operation mode. If False, operation mode is continuous.
        Defaults to `continuous`.
    tau_0 : 
        Cleaning and unloading time (if batch mode). Defaults to 3 hr.
    N :
        Number of reactors.
    heat_exchanger_configuration : 
        What kind of heat exchanger to default to (if any). Valid options include 
        'jacketed', 'recirculation loop', and 'internal coil'. Defaults to 'recirculation loop'.
    dT_hx_loop : 
        Maximum change in temperature for the heat exchanger loop. Defaults to 5 K.
    jacket_annular_diameter :
        Annular diameter of heat exchanger jacket to vessel [m]. Defaults to 0.1 m.
    loading_time :
        Loading time of batch reactor. If not given, it will assume each vessel is constantly
        being filled.
        
    Notes
    -----
    The heat exchanger configuration can be one of the following:

    * 'recirculation loop': 
        The recirculation loop takes into account the required flow rate needed to
        reach the maximum temperature change of the heat exchanger, `dT_hx_loop`. 
        Increasing `dT_hx_loop` decreases the required recirculation flow rate and
        therefore decreases pump costs.
        
        When parallel reactors are required, one recirculation loop (each with a
        pump and heat exchanger) is assumed. Although it is possible to use the
        same recirculation loop for all reactors, this conservative assumption allows
        for each reactor to be operated independently from each other.

    * 'jacketed':
        The jacket does not account for the heat transfer area requirement. 
        It simply assumes that a full jacket can provide the necessary heat transfer
        area to meet the duty requirement. A heuristic annular diameter is assumed
        through `jacket_annular_diameter` (which can be adjusted by the user).
        The temperature at the wall is assumed to be the operating temperature.
        The weight of the jacket is added to the weight of the vessel and the
        cost is compounded together as a jacketed vessel.
        
    * 'internal coil':
        The internal coil is costed as an ordinary helical tube heat exchanger
        with the added assumption that the temperature at the wall is the 
        operating temperature. This method is still not implemented in BioSTEAM
        yet.

    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['H2', 'CO2', 'N2', 'O2', 'H2O', 'AceticAcid'])
    >>> media = bst.Stream(ID='media', H2O=10000, units='kg/hr')
    >>> H2 = bst.Stream(ID='H2', H2=100, units='kg/hr', phase='g')
    >>> fluegas = bst.Stream(ID='fluegas', N2=70, CO2=25, H2O=3, O2=2, units='kg/hr', phase='g')
    >>> # Model acetic acid production from H2 and CO2
    >>> rxn = bst.Rxn('H2 + CO2 -> AceticAcid + H2O', reactant='H2', correct_atomic_balance=True) 
    >>> brxn = rxn.backwards(reactant='AceticAcid')
    >>> R1 = bst.GasFedBioreactor(
    ...     'R1', ins=[media, H2, fluegas], outs=('vent', 'product'), 
    ...     tau=68, V_max=500, controlled_feeds=[1, 2],
    ...     reactions=rxn, backward_reactions=brxn,
    ...     gas_substrates=('H2', 'CO2'),
    ...     titer=dict(AceticAcid=5),
    ...     optimize_power=False,
    ...     kW_per_m3=0.,
    ...     T=305.15
    ... )
    >>> R1.simulate()
    >>> R1.show()
    GasFedBioreactor: R1
    ins...
    [0] media  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  555
    [1] H2  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  17.7
    [2] fluegas  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CO2  1.73
                        N2   7.63
                        O2   0.191
                        H2O  0.509
    outs...
    [0] vent  
        phase: 'g', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2          14.3
                        CO2         0.069
                        N2          7.63
                        O2          0.191
                        H2O         1.09
                        AceticAcid  0.00347
    [1] product  
        phase: 'l', T: 305.15 K, P: 101325 Pa
        flow (kmol/hr): H2          0.00471
                        CO2         0.000672
                        N2          0.00189
                        O2          9.65e-05
                        H2O         556
                        AceticAcid  0.829
    
    """
    _N_ins = 2
    _N_outs = 2
    _ins_size_is_fixed = False
    auxiliary_unit_names = (
        'sparger',
        'compressors',
        'gas_coolers',
        *AbstractStirredTankReactor.auxiliary_unit_names
    )
    T_default = 273.15 + 32 
    P_default = 101325
    kW_per_m3_default = 0.2955 # Reaction in homogeneous liquid
    batch_default = True
    default_methods = AeratedBioreactor.default_methods
    get_kLa = AeratedBioreactor.get_kLa
    get_agitation_power = AeratedBioreactor.get_agitation_power
    
    def _init(self, 
            gas_substrates, 
            reactions=None,
            # Can vary either liquid or gas flows
            titer=None, 
            # Only for controlled gas flows
            backward_reactions=None, 
            controlled_gas_substrates=None,
            controlled_feeds=(), 
            substrate_reactions=None,
            # General design/performance arguments
            design=None, method=None, kLa_kwargs=None,
            theta=0.5, Q_consumption=None,
            cooler_pressure_drop=None,
            compressor_isentropic_efficiency=None,
            **kwargs,
        ):
        if compressor_isentropic_efficiency is None: compressor_isentropic_efficiency = 0.85
        #: Isentropic efficiency of the compressor. Defaults to 0.85.
        self.compressor_isentropic_efficiency = compressor_isentropic_efficiency
        self.cooler_pressure_drop = 20684.28 if cooler_pressure_drop is None else cooler_pressure_drop
        self.theta = theta # Average concentration of gas substrate in the liquid as a fraction of saturation.
        self.Q_consumption = Q_consumption # Forced duty per gas substrate consummed [kJ/kmol].
        self.kLa_kwargs = {} if kLa_kwargs is None else kLa_kwargs
        self.controlled_feeds = controlled_feeds # list[int|Stream] Feed index or stream.
        self.gas_substrates = gas_substrates
        self.titer = titer # dict[str, float] g / L
        self.substrate_reactions = substrate_reactions
        
        controlled_gas_feeds = self.controlled_gas_feeds
        
        if backward_reactions is None and controlled_gas_feeds and titer:
            product = next(iter(titer))
            backward_reactions = (
                reactions.backwards(reactant=product) 
            )
        self.backward_reactions = backward_reactions
        
        if controlled_gas_substrates is None:
            N_gas_feeds = len(controlled_gas_feeds)
            if N_gas_feeds == 0:
                controlled_gas_substrates = ()
            elif N_gas_feeds == len(gas_substrates):
                controlled_gas_substrates = gas_substrates
            else:
                raise ValueError('must specify controlled gas substrates')
        self.controlled_gas_substrates = controlled_gas_substrates
            
        AbstractStirredTankReactor._init(self, reactions=reactions, **kwargs)
        if design is None: 
            if self.kW_per_m3 == 0:
                design = 'Bubble column'
            else:
                design = 'Stirred tank'
        elif design not in aeration.kLa_method_names:
            raise ValueError(
                f"{design!r} is not a valid design; only "
                f"{list(aeration.kLa_method_names)} are valid"
            )
        self.design = design
        if method is None:
            method = self.default_methods[design]
        if (key:=(design, method)) in aeration.kLa_methods:
            self.kLa = aeration.kLa_methods[key]
        elif hasattr(method, '__call__'):
            self.kLa = method
        else:
            raise ValueError(
                f"{method!r} is not a valid kLa method; only "
                f"{aeration.kLa_method_names[design]} are valid"
            )
    
    def _get_duty(self):
        if self.Q_consumption is None:
            H_in = sum(
                [i.H for i in self.ins if i.phase != 'g']
                + [i.outs[0].H for i in self.gas_coolers]
            )
            return self.H_out - H_in + self.Hf_out - self.Hf_in
        else:
            return self.Q_consumption * (
                sum([i.imol['O2'] for i in self.ins])
                - sum([i.imol['O2'] for i in self.outs])
            )
    
    @property
    def vent(self):
        return self._outs[0]
    
    @property
    def controlled_feeds(self):
        return [(i if isinstance(i, bst.Stream) else self.ins[i])
                for i in self._controlled_feeds]
    @controlled_feeds.setter
    def controlled_feeds(self, controlled_feeds):
        self._controlled_feeds = controlled_feeds
    
    @property
    def controlled_liquid_feeds(self):
        return [i for i in self.controlled_feeds if i.phase == 'l']
    
    @property
    def controlled_gas_feeds(self):
        return [i for i in self.controlled_feeds if i.phase == 'g']
    
    @property
    def normal_feeds(self):
        controlled = set(self.controlled_feeds)
        return [i for i in self.ins if i not in controlled]
    
    @property
    def liquid_feeds(self):
        return [i for i in self.ins if i.phase != 'g']
    
    @property
    def gas_feeds(self):
        return [i for i in self.ins if i.phase == 'g']
    
    @property
    def sparged_gas(self):
        return self.sparger-0
    
    def load_auxiliaries(self):
        super().load_auxiliaries()
        self.compressors = []
        self.gas_coolers = []
        for inlet in self.gas_feeds:
            compressor = self.auxiliary(
                'compressors', bst.IsentropicCompressor, inlet, 
                eta=self.compressor_isentropic_efficiency, P=2 * 101325
            )
            self.auxiliary(
                'gas_coolers', bst.HXutility, compressor-0, T=self.T
            )
        self.auxiliary(
            'sparger', bst.Mixer, [i-0 for i in self.gas_coolers]
        )
        
    def get_SURs(self, F_vol):
        produced = bst.Stream(None, thermo=self.thermo)
        for ID, concentration in self.titer.items(): # Titer is in terms of g / 1000 kg Water
            produced.imass[ID] = F_vol * concentration
        consumed = produced.copy()
        self.backward_reactions.force_reaction(consumed)
        SURs = consumed.get_flow('mol/s', self.gas_substrates)
        return SURs
    
    def _get_grouped_substrate_flows(self, streams):
        return self._group_substrate_flows(
            sum([i.get_flow('mol/s', self.gas_substrates) for i in streams])
        ) / 3.6
    
    def _group_substrate_flows(self, F_substrates):
        if self.substrate_reactions is None: return F_substrates
        chemicals = self.chemicals
        flows = chemicals.array(self.gas_substrates, F_substrates)
        self.substrate_reactions.force_reaction(flows)
        return flows[chemicals.get_index(self.controlled_gas_substrates)]
    
    def _run_vent(self, vent, effluent):
        aeration.vent_broth(vent, effluent)
    
    def _run_without_titer_specification(self, effluent, vent, feed, F_substrates=None):
        gas_substrates = self.gas_substrates
        if F_substrates is None: 
            F_substrates = sum([
                i.imol[gas_substrates]
                for i in self.ins
            ]) / 3.6 # mol / s
        
        def f(STR_guess):
            STRs = np.minimum(STR_guess, F_substrates)
            effluent.copy_flow(feed)
            effluent.set_flow(STRs, units='mol/s', key=gas_substrates)
            remaining = F_substrates - STRs
            self._run_reactions(effluent)
            vent.copy_flow(self.sparged_gas)
            vent.imol[gas_substrates] = remaining * 3.6 # mol / s -> kmol / hr
            self._run_vent(vent, effluent) 
            return self.get_STRs()
        
        if np.iscomplex(self.get_STRs()).any():
            breakpoint()
            self.get_STRs()
        self._STRs_last = flx.aitken(f, self.get_STRs(), xtol=1e-6, maxiter=1000)
        
    # def plot_surface_response(self):
    #     self._setup()
        
    #     # Insert code here 
    #     f = lambda x, y: gas_flow_rate_objective(np.array([x, y]))
        
    #     width = 6.6142
    #     aspect_ratio = 0.4
    #     rcParams['figure.figsize'] = (width, width * aspect_ratio)
        
    #     xlim = [guess[0] * 0.2, guess[0] * 5]
    #     ylim = [guess[1] * 0.2, guess[1] * 5]
    #     X, Y, Z = bst.plots.generate_contour_data(
    #         f, xlim=xlim, ylim=ylim, n=20,
    #     )
    #     breakpoint()
    #     # Plot contours
    #     xlabel = "X"
    #     ylabel = 'Y'
    #     metric_bar = bst.plots.MetricBar(
    #         'Error', '', colormaps['viridis_r'],
    #         None, 15, 1
    #     )
    #     fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
    #         X, Y, Z, xlabel, ylabel, None, None, metric_bar,
    #         fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    #     )
            
    def _run(self):
        controlled_gas_feeds = self.controlled_gas_feeds
        controlled_liquid_feeds = self.controlled_liquid_feeds
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        self.sparged_gas.T = vent.T = effluent.T = self.T
        vent.phase = 'g'
        try: feed, = self.liquid_feeds
        except: raise RuntimeError('gas-fed bioreactor must have exactly one liquid feed')
        
        if not self.titer: 
            # Titer is given by the mass transfer; only one substrate is limiting,
            # so we have 1 mass transfer/update rate equation and 1 unknown.
            self._update_liquid_feed()
            self._update_gas_feeds()
            self._run_without_titer_specification(effluent, vent, feed)
            return
        
        if controlled_liquid_feeds and not controlled_gas_feeds:
            # Titer given, must adjust liquid flow so that STR meets SUR.
            # Only one substrate is limiting which gives 1 equation and 1 unknown.
            try: feed, = controlled_liquid_feeds
            except: raise RuntimeError('gas-fed bioreactor must have exactly one liquid feed')
            F_substrates = sum([i.get_flow(units='mol/s', key=self.gas_substrates) for i in self.ins])
            F_liquid_max = self._initialize_controlled_liquid_guess(effluent)
            product, titer = next(iter(self.titer.items()))
            titer /= 1000
            self._update_gas_feeds()
            def liquid_flow_rate_objective(F_feed):
                if F_feed < 0: F_feed = 1e-6
                feed.F_mass = F_feed
                self._update_liquid_feed()
                self._run_without_titer_specification(effluent, vent, feed, F_substrates)
                F_product = vent.imass[product] + effluent.imass[product]
                return F_product / F_feed - titer
            
            flx.aitken_secant(liquid_flow_rate_objective, 0.8 * F_liquid_max, F_liquid_max, ytol=1e-9)
            return
        
        if controlled_gas_feeds and not controlled_liquid_feeds:
            # Controlled gas substrates are the limiting substrates 
            # (this is the # of equations).
            controlled_gas_substrates = self.controlled_gas_substrates
            if len(controlled_gas_feeds) != len(controlled_gas_substrates):
                raise RuntimeError(
                    'number of controlled gas substrates must match controlled '
                    'feeds to close degrees of freedom'
                )
            # Solve gas flow rates to meet titer.
            self._update_liquid_feed()
            effluent.mix_flows(self.liquid_feeds)
            F_vol = feed.F_mass / 1000 # Approximately m3 / hr
            SURs = self._group_substrate_flows(self.get_SURs(F_vol)) # Gas substrate uptake rate [mol / s]
            if (SURs <= 1e-2).all():
                self._run_vent(vent, effluent)
                return
            index = range(len(controlled_gas_substrates))
            baseline_flows = self._get_grouped_substrate_flows(self.normal_feeds)
            substrate_reactions = self.substrate_reactions
            if substrate_reactions is None:
                # Each feed directly controls a gas substrate
                x_substrates = []
                for gas, ID in zip(self.controlled_gas_feeds, controlled_gas_substrates):
                    x_substrates.append(gas.get_molar_fraction(ID))
                
                guess = 1.01 * SURs / x_substrates - baseline_flows
            else:
                # A feed may contribute multiple gas substrates
                controlled_gas_feeds = self.controlled_gas_feeds
                coefficients = np.zeros([len(controlled_gas_substrates), len(controlled_gas_feeds)])
                for i, gas in enumerate(controlled_gas_substrates):
                    for j, stream in enumerate(controlled_gas_feeds):
                        reacted = stream.copy()
                        self.substrate_reactions.force_reaction(reacted)
                        coefficients[i, j] = reacted.imol[gas] / stream.F_mol
                F_min = np.linalg.solve(coefficients, SURs - baseline_flows)
                guess = 1.01 * F_min
            
            def gas_flow_rate_objective(F_feeds):
                F_feeds[F_feeds < 0] *= -1
                for i in index: controlled_gas_feeds[i].set_total_flow(F_feeds[i], 'mol/s')
                self._update_gas_feeds()
                self._run_without_titer_specification(effluent, vent, feed)
                STRs = self._group_substrate_flows(self._STRs_last) # Must meet all substrate demands
                diff = SURs - STRs
                diff[diff < 0] *= 2
                return diff
            
            f = gas_flow_rate_objective
            with catch_warnings():
                filterwarnings('ignore')
                results = fsolve(
                    f, guess, full_output=True, maxfev=500, xtol=1e-9
                )
            self._convergence = results
        elif controlled_liquid_feeds and controlled_gas_feeds: 
            try: feed, = controlled_liquid_feeds
            except: raise RuntimeError('gas-fed bioreactor must have exactly one liquid feed')
            
            # SURs are given by liquid feed and titer (1 equation / 2 unknown)
            # STRs are given by gas feeds (must be N STR/SUR equations and N - 1 unknown feed variables)
            controlled_gas_substrates = self.controlled_gas_substrates
            N_controlled = len(controlled_gas_substrates)
            controlled_feeds = self.controlled_feeds
            if N_controlled != len(controlled_feeds):
                raise RuntimeError(
                    'number of controlled gas substrates must be equal to the number of controlled '
                    'feeds'
                ) # Given there is only one controlled liquid feed, this statement holds true
            
            effluent.copy_flow(feed)
            F_liquid_max = self._initialize_controlled_liquid_guess(effluent, maxflow=True)
            feed.F_mass = F_liquid_max
            SURs = self._group_substrate_flows(self.get_SURs(F_liquid_max / 1000)) # Gas substrate uptake rate [mol / s]
            substrate_reactions = self.substrate_reactions
            baseline_flows = self._get_grouped_substrate_flows(self.normal_feeds)
            if substrate_reactions is None:
                # Each feed directly controls a gas substrate
                x_substrates = []
                subset = []
                for i, (stream, gas) in enumerate(zip(controlled_feeds, controlled_gas_substrates)):
                    if stream.phase != 'g': continue
                    subset.append(i)
                    x_substrates.append(gas.get_molar_fraction(ID))
                F_min = SURs[subset] / x_substrates - baseline_flows[subset]
            else:
                # A feed may contribute multiple gas substrates
                coefficients = []
                subset = []
                for gas in controlled_gas_substrates:
                    row = []
                    coefficients.append(row)
                    for stream in controlled_feeds:
                        if stream.phase != 'g': continue
                        reacted = stream.copy()
                        self.substrate_reactions.force_reaction(reacted)
                        row.append(reacted.imol[gas] / stream.F_mol)
                for i, (stream, gas) in enumerate(zip(controlled_feeds, controlled_gas_substrates)):
                    if stream.phase != 'g': continue
                    subset.append(i)
                coefficients = np.array(coefficients)
                F_min = np.linalg.solve(coefficients[subset], SURs[subset] - baseline_flows[subset]) 
                
            index = range(N_controlled - 1)
            def gas_flow_rate_objective(F_controlled):
                F_controlled[F_controlled < 0] *= -1 
                for i in index:
                    gas = controlled_gas_feeds[i]
                    gas.set_total_flow(F_controlled[i], 'mol/s')
                feed.F_mass = F_liq = F_controlled[-1]
                F_substrates = sum([
                    i.imol[self.gas_substrates]
                    for i in self.ins
                ]) / 3.6 # mol / s
                self._update_liquid_feed()
                self._update_gas_feeds()
                self._run_without_titer_specification(effluent, vent, feed, F_substrates)
                SURs = self._group_substrate_flows(self.get_SURs(F_liq / 1000)) # Gas substrate uptake rate [mol / s]
                STRs = self._group_substrate_flows(self._STRs_last) # Must meet all substrate demands
                return SURs - STRs
            
            f = gas_flow_rate_objective
            guess = np.array([*F_min, F_liquid_max])
            with catch_warnings():
                filterwarnings('ignore')
                results = fsolve(
                    f, guess, full_output=True, maxfev=500, xtol=1e-9
                )
            self._convergence = results
        else:
            raise RuntimeError('cannot satisfy titer specification without controlled feeds')
        
    def _run_reactions(self, effluent, maxflow=None):
        if self.substrate_reactions: self.substrate_reactions.force_reaction(effluent)
        if maxflow: 
            data = effluent.get_data()
            rxns = self.reactions
            rxns.force_reaction(effluent)
            for reactant in self.gas_substrates:
                if effluent.imol[reactant] > 0:
                    if isinstance(rxns, bst.Rxn):
                        rxns.reactant = reactant
                    else:
                        for rxn in rxns:
                            if rxn.istoichiometry[reactant] < 0:
                                rxn.reactant = reactant
                
                    effluent.set_data(data)
                    rxns.force_reaction(effluent)
                    break
        else:
            data = effluent.get_data()
            rxns = self.reactions
            rxns.force_reaction(effluent)
            for reactant in self.gas_substrates:
                if effluent.imol[reactant] < 0:
                    if isinstance(rxns, bst.Rxn):
                        rxns.reactant = reactant
                    else:
                        for rxn in rxns:
                            if rxn.istoichiometry[reactant] < 0:
                                rxn.reactant = reactant
                
                    effluent.set_data(data)
                    rxns.force_reaction(effluent)
                    break
        
    def _initialize_controlled_liquid_guess(self, effluent, maxflow=None):
        effluent.mix_flows(self.ins)
        self._run_reactions(effluent, maxflow)
        product, titer = next(iter(self.titer.items()))
        F_liquid_max = 1000 * effluent.imass[product] / titer # kg / hr
        return F_liquid_max
        
    def _solve_total_power(self, SURs): # For STR = SUR [mol / s]
        gas_in = self.sparged_gas
        N_reactors = self.parallel['self']
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        D = self.get_design_result('Diameter', 'm')
        F = gas_in.get_total_flow('m3/s') / N_reactors / operating_time
        R = 0.5 * D
        A = pi * R * R
        self.superficial_gas_flow = U = F / A # m / s 
        vent = self.vent
        Ps = []
        for gas_substrate, SUR in zip(self.gas_substrates, SURs):
            Py_gas = gas_in.get_property('P', 'bar') * gas_in.imol[gas_substrate] / gas_in.F_mol
            Py_vent = 0. if vent.isempty() else vent.get_property('P', 'bar') * vent.imol[gas_substrate] / vent.F_mol
            C_sat_gas = aeration.C_L(self.T, Py_gas, gas_substrate) # mol / kg
            C_sat_vent = aeration.C_L(self.T, Py_vent, gas_substrate) # mol / kg
            theta = self.theta
            LMDF = aeration.log_mean_driving_force(C_sat_vent, C_sat_gas, theta * C_sat_vent, theta * C_sat_gas)
            kLa = SUR / (LMDF * V * self.effluent_density * N_reactors * operating_time)
            Ps.append(aeration.P_at_kLa_Riet(kLa, V, U, **self.kLa_kwargs))
        P = max(Ps)  
        agitation_power_kW = P / 1000
        compressor_power_kW = sum([i.power_utility.consumption for i in self.compressors]) / N_reactors
        total_power_kW = (agitation_power_kW + compressor_power_kW) / V
        self.kW_per_m3 = agitation_power_kW / V 
        return total_power_kW
    
    def get_STRs(self):
        """Return the gas substrate transfer rate in mol/s."""
        V = self.get_design_result('Reactor volume', 'm3') * self.V_wf
        operating_time = self.tau / self.design_results.get('Batch time', 1.)
        N_reactors = self.parallel['self']
        gas_in = self.sparged_gas
        kLa = self.get_kLa() # 1 / s 
        vent = self.vent
        P_gas = gas_in.get_property('P', 'bar')
        P_vent = vent.get_property('P', 'bar')
        STRs = []
        for ID in self.gas_substrates:
            Py_gas = P_gas * gas_in.imol[ID] / gas_in.F_mol
            Py_vent = P_vent * vent.imol[ID] / (vent.F_mol or 1)
            C_sat_gas = aeration.C_L(self.T, Py_gas, ID) # mol / kg
            C_sat_vent = aeration.C_L(self.T, Py_vent, ID) # mol / kg
            theta = self.theta
            LMDF = aeration.log_mean_driving_force(C_sat_vent, C_sat_gas, theta * C_sat_vent, theta * C_sat_gas)
            STRs.append(
                kLa * LMDF * self.effluent_density * V * N_reactors * operating_time # mol / s
            )
        return np.array(STRs)
    
    def _update_liquid_feed(self):
        AbstractStirredTankReactor._design(self, size_only=True)
        liquid = bst.Stream(None, thermo=self.thermo)
        liquid.mix_flows([i for i in self.ins if i.phase != 'g'])
        liquid.copy_thermal_condition(self.outs[0])
        self.effluent_density = rho = liquid.rho
        length = self.get_design_result('Length', 'm') * self.V_wf
        P_inlet = g * rho * length + 101325
        P_compressor = P_inlet + self.cooler_pressure_drop # Pa
        for i in self.compressors: i.P = P_compressor
        self.sparged_gas.P = P_inlet
        
    def _update_gas_feeds(self):
        self.sparged_gas.mix_flows(self.gas_feeds)
        
    def _design(self):
        for i in self.compressors: i.simulate()
        for i in self.gas_coolers: i.simulate()
        self.sparger.simulate()
        AbstractStirredTankReactor._design(self)
        self.parallel['sparger'] = 1
        self.parallel['compressors'] = 1
        self.parallel['gas_coolers'] = 1
        # For robust process control, do not include in HXN
        for unit in self.auxiliary_units:
            for hu in unit.heat_utilities: hu.hxn_ok = False
    
GFB = GasFedBioreactor