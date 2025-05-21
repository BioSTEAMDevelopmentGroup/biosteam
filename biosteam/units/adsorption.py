# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. contents:: :local:

.. autoclass:: biosteam.units.adsorption.SingleComponentAdsorptionColumn

References
----------
.. [1] Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE

.. [2] Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).

.. [3] Seider, W. D., Lewin,  D. R., Seader, J. D., Widagdo, S., Gani,
    R., & Ng, M. K. (2017). Product and Process Design Principles. Wiley.
    Cost Accounting and Capital Cost Estimation (Chapter 16)

"""
import biosteam as bst
from .design_tools import PressureVessel
from numpy.testing import assert_allclose
import numpy as np
from numba import njit
import flexsolve as flx
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.ndimage.filters import gaussian_filter
from thermosteam.units_of_measure import format_units

__all__ = ('SingleComponentAdsorptionColumn', 'AdsorptionColumn',)

@njit(cache=True)
def fixed_step_odeint(f, x0, step, tf, args):
    M = int(tf/step)
    N = len(x0)
    ts = np.linspace(0, tf, M)
    xs = np.zeros((M, N))
    xs[0] = x0
    for i in range(1, M):
        dx_dt = f(ts[i], x0, *args)
        x0 += dx_dt * step
        xs[i] = x0 
    return (ts, xs)

@njit(cache=True)
def equilibrium_loading_Langmuir_isotherm(
        C,  # Solute concentration in the feed [g / L]
        KL,  # Equilibrium constant
        q_max,  # Maximum equilibrium loading mg/g
    ):
    return C * KL * q_max / (1 + KL * C)  # g / kg

@njit(cache=True)
def equilibrium_loading_Freundlich_isotherm(
        C,  # Solute concentration in the feed [g / L]
        KL,  
        n,  
    ):
    return KL * C**(1/n) # g / kg

@njit(cache=True)
def estimate_equilibrium_bed_length(
        C,  # Solute concentration in the feed [kg / m3]
        u,  # Superficial velocity [m / h]
        cycle_time,  # Cycle time [h]
        rho_adsorbent,  # Bulk density of the bed [kg / m3]
        q0,  # Equilibrium loading [g / kg]
    ):
    q_capacity = q0 * rho_adsorbent / 1000 # kg / m3
    L = C * u * cycle_time / q_capacity  # m
    return L

@njit(cache=True)
def dCdt_optimized_Langmuir(
        t, C, N_slices, 
        Da, dz, KL_qm,
        q0_over_C0,
        q0_KL,
        beta_q0_rho_over_C0,
    ):
    CL = C[:N_slices] 
    qt = C[N_slices:] # q-avg
    CL_ = CL.copy()
    CL_[CL < 0] = 0
    qe = (CL_ * KL_qm)/(q0_over_C0 + q0_KL * CL_)
    dCL_dz = CL_
    dCL_dz[1:] = (CL[1:] - CL[:-1]) / dz
    dCL_dz[0] = (CL[0] - 1) / dz
    dC_dt = C.copy()
    dq_dt = Da * (qe - qt)
    dCL_dt = -dCL_dz - beta_q0_rho_over_C0 * dq_dt
    dC_dt[:N_slices] = dCL_dt
    dC_dt[N_slices:] = dq_dt
    return dC_dt 

@njit(cache=True)
def dCdt(
        t, C, N_slices, 
        Da, dz, isotherm_model, isotherm_args,
        beta_q0_rho_over_C0, C0, q0
    ):
    CL = C[:N_slices] 
    qt = C[N_slices:] # q-avg
    CL_ = CL * C0
    mask = CL_ < 0
    CL_[mask] = 0
    qe = isotherm_model(CL_, *isotherm_args) / q0
    qe[qe > q0] = q0
    dCL_dz = CL.copy()
    dCL_dz[1:] = (CL[1:] - CL[:-1]) / dz
    dCL_dz[0] = (CL[0] - 1) / dz
    dC_dt = C.copy()
    dq_dt = Da * (qe - qt)
    dCL_dt = -dCL_dz - beta_q0_rho_over_C0 * dq_dt
    dC_dt[:N_slices] = dCL_dt
    dC_dt[N_slices:] = dq_dt
    return dC_dt 

def adsorption_bed_pressure_drop(
        D = 0.001, # particle diameter [m]
        rho = 800, # liquid density [kg/m3]
        mu = .001, # viscosity [kg/m*s]
        epsilon = 1 - 0.38/0.8, # void fraction
        u = 14/3600, # superficial velosity [m/s]
        L = 5, # length of bed [m]
    ):
    re = 1 - epsilon
    e2 = epsilon * epsilon
    e3 = e2 * epsilon
    N_Re = D * u * rho/(mu*re)
    if N_Re < 10: # Kozeny-Carman equation
        dP = 150 * mu * u * re ** 2 * L/(D ** 2 * e3)
    elif 10 < N_Re < 1000: # Ergun equation 
        dP = 1.75 * u**2 * L * re * rho/(D * e3) + 150 * mu * u * re**2/(D**2 * e3)
    else: # Burke-Plummer equation
        dP = 1.75 * u**2 * L * re * rho/(D * e3) + 150
    return dP

class SingleComponentAdsorptionColumn(PressureVessel, bst.Unit):
    """
    Create an adsorption column which can operate by temperature or pressure swing,
    or by disposing the adsorbent. Heats of adsorption are neglected but may become
    an option with future development in BioSTEAM. Default parameters are heuristic values 
    for adsorption of water and polar components using activated carbon. The
    design and cost algorithm of the adsorption column is based on [1]_, [2]_, 
    and [3]_. 
    
    By default, three vessels are used: a lead, a guard, and a standby vessel. 
    The lead is in equilibrium, guard holds the breakthrough curve and mass 
    transfer zone, and the standby is in regeneration. Once the standby vessel
    is regenerated, it becomes the guard, the lead becomes the standby, and 
    the guard becomes the lead. Therefore, the cycle time is the regeneration
    time. 
    
    For temperature or pressure swing adsorption, set the temperature or pressure
    of the regeneration fluid. If a regeneration fluid is used, the flow rate 
    of regeneration fluid solved for such that it meets the required cycle time.
    
    Parameters
    ----------
    ins : 
        * [0] Feed
        * [1] Regeneration fluid
        * [2] Adsorbent (if no regeneration fluid is used).
    outs : 
        * [0] Effluent
        * [1] Purge
        * [2] Spent adsorbent (if no regeneration fluid is used)
    k : float, optional
        Mass transfer coefficient [1/h]. If not given, mass transfer model is 
        ignored.
    superficial_velocity : float, optional
        Superficial velocity of the feed. The diameter of the receiving vessel adjusts
        accordingly. Defaults to 7.2 [m / hr]. Typical velocities are 4 to 14.4 m / hr for liquids [1]_.
        Common velocity range for gasses is 504-2160 m / hr [1]_.
    cycle_time : float, optional
        Time at which the receiving vessel is switched. Defaults to 3 [hr]. 
        Note that 1-2 hours required for thermal-swing-adsorption 
        (TSA) for silica gels [2]_. 
    isotherm_args : tuple[float, ...], optional
        Arguments to pass to the isotherm model. 
        If isotherm model is 'Langmuir', arguments must be:
        * KL: Equilibrium constant [L/mg] for Langmuir isotherm
        * q_max: Maximum equilibrium loading [mg/g] for Langmuir isotherm
        If isotherm model is 'Freundlich', arguments must be:
        * K: Equilibrium constant [L/mg] for Freundlich isotherm
        * n: exponential coefficient for Freundlich isotherm
    isotherm_model : Callable|str, optional, 
        Can be 'Langmuir', 'Freundlich', or a function.
        If a function is given, it should be in the form of f(C, *args) -> g absorbate / kg absorbent.
    void_fraction : 0.525 
        Fraction of empty space in the adsorbent by vol.
    rho_adsorbent : float, optional
        The density of the adsorbent. Defaults to 380 [kg/m3], which is common 
        for activated carbon granules.
    regeneration_fluid : dict[str, float]
        Arguments to initialize fluid used to regenerate the bed.
    LUB : float, optional 
        Length of unused bed. Defaults to 1.219 m if no mass transfer 
        coefficient is given. This parameter is ignored if isotherms are provided.
    P : float, optional
        Operating pressure of the column. Defaults to 101325.
    C_final_scaled : float, optional
        Final outlet concentration at breakthrough relative to inlet. 
        Defaults to 0.05.
    regeneration_fluid : dict[str, float]
        Regeneration fluid composition and thermal conditions.
    k_regeneration : float, optional 
        Mass transfer coefficient of the regeration solvent [1/h].
    regeneration_isotherm_args : tuple[float, ...], optional
        Arguments to regeneration isotherm model.
    regeneration_isotherm_model : Callable|str
        Isotherm model for regenerating the adsorbent.
    N_slices : int, optional
        Number of slices to model mass transfer. Defaults to 80.
    adsorbate : str
        Name of adsorbate.
    vessel_material : float, optional
        Vessel material. Defaults to 'Stainless steel 316',
    vessel_type : float, optional
        Vessel type. Defaults to 'Vertical'.
    adsorbent : str
        Name of adsorbent. Defaults to 'ActivatedCarbon'.
    N_columns : int
        Number of columns can be either 2 or 3. The 3-column configuration 
        minimizes cost of adsorbent regeneration. The 2-column configuration 
        minimizes capital cost.
    particle_diameter : float
        Particle diameter assumed for pressure drop estimation [m]. 
        Defaults to 0.004 m.    
    
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo([
    ...    'Water', 
    ...    bst.Chemical('Adsorbate', search_db=False, default=True, phase='l'),
    ...    bst.Chemical('ActivatedCarbon', search_db=False, default=True, phase='s')
    ... ])
    >>> feed = bst.Stream(ID='feed', phase='l', T=298, P=1.01e+06,
    ...                   Water=1000, Adsorbate=0.001, units='kg/hr')
    >>> A1 = bst.AdsorptionColumn(
    ...     'A1', ins=feed, outs='effluent',
    ...     cycle_time=1000,
    ...     superficial_velocity=9.2,
    ...     isotherm_model='Langmuir',
    ...     isotherm_args=(1e3, 7.), # K [g / L], q_max [g / kg]
    ...     k=0.3, # [1 / hr]
    ...     void_fraction=0.525, 
    ...     C_final_scaled=0.05,
    ...     adsorbate='Adsorbate',
    ...     particle_diameter=0.004,
    ... )
    >>> A1.simulate()
    >>> # A1.simulation_gif() # Create a gif of the fluid concentration profile over time.
    
    .. image:: ../../images/adsorption_column.gif
       :align: center
       :alt: Adsorption column fluid concentration profile over time.
    
    >>> A1.results()
    Single component adsorption column                  Units                   A1
    Electricity         Power                              kW                0.241
                        Cost                           USD/hr               0.0188
    Design              Diameter                           ft                 1.22
                        Length                             ft                 12.2
                        Vessel type                                       Vertical
                        Weight                             lb                  521
                        Wall thickness                     in                 0.25
                        Pressure drop                      Pa                  103
                        Vessel material                        Stainless steel 316
    Purchase cost       Vertical pressure vessel (x3)     USD             6.12e+04
                        Platform and ladders (x3)         USD             8.37e+03
                        Pump - Pump (x3)                  USD             1.31e+04
                        Pump - Motor (x3)                 USD                  174
    Total purchase cost                                   USD             8.28e+04
    Utility cost                                       USD/hr               0.0188
    
    >>> A1.show('wt')
    SingleComponentAdsorptionColumn: A1
    ins...
    [0] feed  
        phase: 'l', T: 298 K, P: 1.01e+06 Pa
        flow (kg/hr): Water      1e+03
                      Adsorbate  0.001
    [1] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): ActivatedCarbon  0.154
    outs...
    [0] effluent  
        phase: 'l', T: 298 K, P: 101325 Pa
        flow (kg/hr): Water  1e+03
    [1] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kg/hr): ActivatedCarbon  0.154
    
    """
    auxiliary_unit_names = (
        'pump', 'heat_exchanger',
        'regeneration_pump'
    )
    _N_ins = 3
    _N_outs = 3

    _units = {'Pressure drop': 'Pa',
              **PressureVessel._units}

    # Cost of regeneration in $/m3
    adsorbent_cost = {
        'ActivatedAlumina': 72 * 0.0283168,
        'ActivatedCarbon': 41 * 0.0283168,
        'SilicaGel': 210 * 0.0283168,
        'MolecularSieves': 85 * 0.0283168,
    }

    # TODO: Update later with reference
    # This is only active if modeling regeneration of solute.
    _default_equipment_lifetime = {
        'ActivatedAlumina': 100,
        'ActivatedCarbon': 100,
        'SilicaGel': 100,
        'MolecularSieves': 100,
    }

    isotherm_models = {
        'Langmuir': equilibrium_loading_Langmuir_isotherm,
        'Freundlich': equilibrium_loading_Freundlich_isotherm,
    }

    def _init(self,
            cycle_time,
            # If not given, mass transfer model is ignored
            k=None,  # Mass transfer coefficient [1/h]
            # Langmuir:
            # - KL: Equilibrium constant [L/mg] for Langmuir isotherm
            # - q_max: Maximum equilibrium loading [mg/g] for Langmuir isotherm
            # Freundlich:
            # - K: Equilibrium constant [L/mg] for Freundlich isotherm
            # - n: exponential coefficient for Freundlich isotherm
            isotherm_args=None,
            # Can be a function, 'Langmuir', or 'Freundlich'.
            # If a function is given, it should be in the form of f(C, *args) -> g absorbate / kg absorbent
            isotherm_model=None, 
            k_regeneration=None, # Mass transfer coefficient of the regeration solvent [1/h]
            regeneration_isotherm_args=None,
            regeneration_isotherm_model=None,
            # Note that the columns are sized according to the limiting isotherm.
            LUB_forced=None, # Defaults to 0.6096 m if no mass transfer coefficient is given.
            void_fraction=1 - 0.38 / 0.8, # Solid by vol [%]
            rho_adsorbent=380, # kg/m3
            superficial_velocity=14.4,  # [m / h]
            P=101325,
            C_final_scaled=0.05, # Final outlet concentration at breakthrough relative to inlet.
            N_slices=int(50), # Number of slices to model mass transfer.
            regeneration_fluid=None, # [dict] Regeneration fluid composition and thermal conditions.
            adsorbate=None, # [str] Name of adsorbate.
            vessel_material='Stainless steel 316',
            vessel_type='Vertical',
            adsorbent='ActivatedCarbon',
            N_columns=3, # 3 column configuration minimizes cost of adsorbent regeneration. 2 column configuration minimizes capital cost.
            particle_diameter = 0.004, # [m]
        ):
        if N_columns not in (2, 3):
            raise ValueError('only 2 or 3-column configurations are valid')
        if isotherm_model is None:
            if isotherm_args is not None:
                raise ValueError('no isotherm model given')
            isotherm_model = lambda C: 0.1
        elif isinstance(isotherm_model, str):
            if isotherm_model in self.isotherm_models:
                isotherm_model = self.isotherm_models[isotherm_model]
            else:
                raise ValueError(
                    f'isotherm model {isotherm_model} is not available; '
                    f'valid options are {list(self.isotherm_models)}'
                )
        if regeneration_isotherm_model is None:
            if regeneration_isotherm_args is not None:
                raise ValueError('no regeneration isotherm model given')
        elif isinstance(isotherm_model, str):
            if regeneration_isotherm_model in self.isotherm_models:
                regeneration_isotherm_model = self.isotherm_models[regeneration_isotherm_model]
            else:
                raise ValueError(
                    f'isotherm model {regeneration_isotherm_model} is not available; '
                    f'valid options are {list(self.isotherm_models)}'
                )
        if k is None and k_regeneration is None:
            if LUB_forced is None: LUB_forced = 0.6096
        elif LUB_forced is not None:
            raise ValueError('length of unused bed given, but will be '
                             'estimated based on mass transfer modeling')
        if regeneration_isotherm_args is None: regeneration_isotherm_args = ()
        if isotherm_args is None: isotherm_args = ()
        self.adsorbate = adsorbate
        self.regeneration_fluid = regeneration_fluid
        self.N_columns = N_columns
        self.cycle_time = cycle_time
        self.superficial_velocity = superficial_velocity
        self.P = P
        self.isotherm_args = isotherm_args
        self.isotherm_model = isotherm_model
        self.regeneration_isotherm_args = regeneration_isotherm_args
        self.regeneration_isotherm_model = regeneration_isotherm_model
        self.k = k
        self.k_regeneration = k_regeneration
        self.void_fraction = void_fraction
        self.rho_adsorbent = rho_adsorbent
        self.C_final_scaled = C_final_scaled
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.N_slices = N_slices
        self.particle_diameter = particle_diameter
        self.adsorbent = adsorbent
        self.LUB_forced = LUB_forced
        self.auxiliary('pump', bst.Pump, ins=self.ins[0])
        if regeneration_fluid:
            regeneration_pump = self.auxiliary('regeneration_pump', bst.Pump, ins=self.ins[1])
            self.auxiliary('heat_exchanger', bst.HXutility, ins=regeneration_pump.outs[0])

    def _run(self):
        feed, regeneration_fluid, adsorbent = self.ins
        outlet, spent_fluid, spent_adsorbent = self.outs
        regeneration_fluid.empty()
        spent_fluid.empty()
        adsorbent.empty()
        spent_adsorbent.empty()
        outlet.copy_like(feed)
        outlet.P = self.P
        outlet.imol[self.adsorbate] = 0
        if self.regeneration_fluid:
            regeneration_fluid = self.heat_exchanger.outs[0]
            regeneration_fluid.reset(**self.regeneration_fluid)
            self.heat_exchanger.T = regeneration_fluid.T
            self._size_columns()
            self.regeneration_superficial_velocity = u = self._solve_regeneration_velocity()
            Q = self.area * u
            regeneration_fluid.F_vol = Q * 3600 # m3/s to m3/h
            self.heat_exchanger.ins[0].copy_flow(regeneration_fluid)
            self.heat_exchanger.run()
            spent_fluid.copy(regeneration_fluid)
            spent_fluid.imol[self.adsorbate] = feed.imol[self.adsorbate]

    def column_widget(self):
        # Not currently working
        import ipywidgets as widgets
        from ipywidgets import interact
        t = self.t
        N_time = len(t)
        tf = t[-1]
        z = np.arange(0, 1, 1/self.N_slices)
        CL_scaled = self.CL_scaled
        plt.xlabel('z')
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.ylabel('Concentration')
        def plot(time):
            time_index = int(time / tf * N_time) - 1
            plt.clf()
            plt.plot(z, CL_scaled[:, time_index], label = '$\theta$ vs z')
            plt.show()
        
        slider = widgets.FloatSlider(min=0, max=tf, step=tf/100, value=0.5*tf)
        return interact(plot, time=slider)

    def simulation_gif(self, liquid=True):
        import matplotlib.animation as animation
        from scipy.interpolate import interp1d
        t = self.t_scaled
        z = np.arange(0, 1, 1/self.N_slices)
        y = self.CL_scaled if liquid else self.q_scaled
        fig = plt.figure()
        N_frames = 120
        f = interp1d(t, y)
        x = np.linspace(0, t[-1], N_frames)
        y = f(x)
        def animate(frame):
            plt.clf()
            artist = plt.plot(z, y[:, frame])
            plt.xlabel(r'z = x / L')
            plt.ylabel(r'$\theta$ = C / C$_{feed}$')
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
            plt.ylim([0, 1])
            return (artist,)
        
        ani = animation.FuncAnimation(
            fig, animate, repeat=True,
            frames=N_frames - 1, interval=66,
        )
        
        # To save the animation using Pillow as a gif
        writer = animation.PillowWriter(fps=15,
                                        metadata=dict(artist='Me'),
                                        bitrate=1800)
        ani.save('adsorption_column.gif', writer=writer)
        return ani

    def _simulate_adsorption_bed(
            self,
            regeneration,
            L, # m
            u, # [m / hr]
        ):
        cycle_time = self.cycle_time
        N_slices = self.N_slices
        C_init = np.zeros(2 * N_slices)
        void_fraction = self.void_fraction
        rho_adsorbent = self.rho_adsorbent
        if regeneration:
            C_init[N_slices:] = 1 # saturated column for regeneration
            isotherm_model = self.regeneration_isotherm_model
            isotherm_args = self.regeneration_isotherm_args
            k = self.k_regeneration
        else:
            k = self.k
            isotherm_model = self.isotherm_model
            isotherm_args = self.isotherm_args
        Da = k * L / u
        beta = (1 - void_fraction)/void_fraction
        self.dz = dz = 1 / (N_slices - 1)
        t_scale = (L / u)
        C0 = self.C0
        q0 = isotherm_model(C0, *isotherm_args)
        beta_q0_rho_over_C0 = (beta * q0 * rho_adsorbent / 1000) / C0
        if isotherm_model is equilibrium_loading_Langmuir_isotherm:
            KL, q_m = isotherm_args
            KL_qm = KL * q_m
            q0_over_C0 = q0 / C0
            q0_KL = q0 * KL
            args = (N_slices, Da, dz, KL_qm,
                    q0_over_C0, q0_KL, beta_q0_rho_over_C0)
            f = dCdt_optimized_Langmuir
            sol = solve_ivp(
                f, t_span=(0, cycle_time / t_scale), 
                y0=C_init, 
                method='BDF',
                args=args,
            ) 
            CL = sol.y[:N_slices]
            q = sol.y[N_slices:]
            t = sol.t
        else:
            args = (N_slices, Da, dz, isotherm_model, isotherm_args,
                    beta_q0_rho_over_C0, C0, q0)
            f = dCdt
            time_step = 0.2 / N_slices
            t, Y = fixed_step_odeint(
                f, C_init, time_step, cycle_time / t_scale, args
            )
            t = np.array(t)
            Y = gaussian_filter(np.array(Y).T, 5, axes=1)
            CL = Y[:N_slices]
            q = Y[N_slices:]
        if regeneration:
            self.rt_scaled = t
            self.rCL_scaled = CL
            self.rq_scaled = q
            self.rt_scale = t_scale
        else:
            self.t_scaled = t
            self.CL_scaled = CL
            self.q_scaled = q
            self.t_scale = t_scale
    
    def _estimate_length_of_unused_bed(self, LUB):
        L_guess = self.LES + LUB
        self._simulate_adsorption_bed(
            False, L_guess, self.superficial_velocity, 
        )
        CL = self.CL_scaled
        C_min = self.C_final_scaled
        C_max = 1 - C_min
        mask = (CL > C_max).any(axis=0) & (CL < C_min).any(axis=0) # slices that have breakthrough curve
        time_index = np.where(mask)[0][-1] # As close to breakthrough as possible
        start_index = np.where(CL[:, time_index] >= C_max)[0][-1]
        end_index = np.where(CL[:, time_index] <= C_min)[0][0]
        length = end_index - start_index
        self.t_scaled = self.t_scaled[:time_index]
        self.CL_scaled = self.CL_scaled[:, :time_index]
        self.q_scaled = self.q_scaled[:, :time_index]
        self.MTZ = length * self.dz * L_guess
        self.LUB = LUB = self.MTZ / 2
        return LUB
    
    def _regeneration_time_objective(self, 
            u, # [m / hr]
        ):
        self._simulate_adsorption_bed(
            True, self.column_length, u,
        )
        q = self.q_scaled
        mask = (q > 1e-6).all(axis=0) & (q < 1e-3).all(axis=0) # slices that have finished regenerating
        time_index = np.where(mask)[0][-1] # Conservative time point
        return self.t_scaled[time_index] * self.t_scale - self.cycle_time
    
    def _solve_regeneration_velocity(self):
        self.regeneration_velocity = u = flx.IQ_interpolation(
            self._estimate_regeneration_time, 4, 14.4, self.cycle_time, 
            xtol=1e-3, ytol=1e-2, 
        ) 
        return u
        
    # def _final_concentration_objective(self, 
    #         L_guess, # m
    #         u, # [m / hr]
    #         cycle_time, # [h]
    #         void_fraction, # external void fraction
    #         q_m, # [g / kg]
    #         k, # [1 / hr]
    #         KL, # [L / g]
    #         feed_rho, # [kg / L]     
    #         C0, # [g / L]
    #     ):
    #     self._simulate_adsorption_bed(
    #         L_guess, u, cycle_time, void_fraction, q_m, k, KL, feed_rho, C0, 
    #     )
    #     C_final = self.CL_scaled[-1, -1] # Final concentration at the end of the bed after cycle time
    #     return C_final - self.C_final
    
    # def _solve_length_of_unused_bed(self, 
    #         u, # [m / hr]
    #         cycle_time, # [h]
    #         void_fraction, # external void fraction
    #         q_m, # [g / kg]
    #         k, # [1 / hr]
    #         KL, # [L / g]
    #         feed_rho, # [kg / L]     
    #         C0, # [g / L]
    #     ):
        
    #     L_guess = self.LES + getattr(self, 'LUB', 2/3) # # 2/3 is a heuristic from AICHE adsorption basics part 1
    #     L = secant(
    #         self._final_concentration_objective, L_guess, 
    #         xtol=self.LES/100, ytol=self.C_final/100,
    #         args=(u, cycle_time, void_fraction, q_m, k, KL, feed_rho, C0)
    #     ) 
    #     LUB = L - self.LES
    #     assert LUB > 0
    #     return LUB

    def estimate_length_of_unused_bed(self):
        LUB_guess = getattr(self, 'LUB', 2/3) # 2/3 is a heuristic from AICHE adsorption basics part 1
        return flx.wegstein(
            self._estimate_length_of_unused_bed, LUB_guess,
            xtol=self.LES * 1e-2, maxiter=3, checkiter=False,
        )

    def _size_columns(self):
        # Design based on capacity (as estimated by the Langmuir isotherm).
        # The length of the mass transfer zone (MTZ) is assumed to be 4 ft based
        # on AICHE adsorption basics part 1. However, future work will use the rate
        # constant and mass transfer modeling to estimate the length of the MTZ.
        feed = self.ins[0]
        Q = feed.F_vol  # m3 / hr
        u = self.superficial_velocity
        self.area = A = Q / u
        self.diameter = diameter = 2 * (A/np.pi) ** 0.5
        self.C0 = C_feed = feed.imass[self.adsorbate] / Q  # g/L or kg/m3
        if C_feed == 0: 
            self.q0 = self.LES = self.LUB = self.total_length = self.column_length = 0
            return 0, 0
        cycle_time = self.cycle_time
        self.q0 = q0 = self.isotherm_model(
            C_feed, *self.isotherm_args
        )
        self.LES = LES = estimate_equilibrium_bed_length(
            C_feed, u, cycle_time, self.rho_adsorbent, q0
        )
        if self.LUB_forced is None:
            self.LUB = LUB = self.estimate_length_of_unused_bed()
        else:
            self.LUB = LUB = self.LUB_forced
        self.total_length = total_length = LES + LUB
        if self.N_columns == 3: 
            column_length = total_length / 2
        elif self.N_columns == 2:
            column_length = total_length
        self.column_length = column_length
        return diameter, column_length

    def _design(self):
        diameter, column_length = self._size_columns()
        feed = self.ins[0]
        if column_length == 0:
            self.design_results['Weight'] = 0
            self.design_results['Diameter'] = 0
            self.design_results['length'] = 0
            return
        self.set_design_result('Diameter', 'm', diameter)
        self.set_design_result('Length', 'm', column_length)
        self.design_results.update(
            self._vertical_vessel_design(
                feed.get_property('P', 'psi'),
                self.design_results['Diameter'],
                self.design_results['Length']
            )
        )
        dP = adsorption_bed_pressure_drop(
            D = self.particle_diameter,
            rho = feed.rho,
            mu = feed.get_property('mu', 'kg/m/s'), # viscosity [kg/m*s]
            epsilon = self.void_fraction, # void fraction
            u = self.superficial_velocity / 3600, # superficial velosity [m/s]
            L = self.total_length, # length of bed [m]
        )
        self.set_design_result('Pressure drop', 'Pa', dP)
        self.pump.P = (self.P - feed.P) + dP
        self.pump.simulate()
        if self.regeneration_fluid:
            feed = self.ins[1]
            outlet = self.outs[1]
            self.regeneration_pump.P = (outlet.P - feed.P) + adsorption_bed_pressure_drop(
                D = self.particle_diameter,
                rho = feed.rho,
                mu = feed.get_property('mu', 'kg/m/s'), # viscosity [kg/m*s]
                epsilon = self.void_fraction, # void fraction
                u = self.regeneration_superficial_velocity / 3600, # superficial velosity [m/s]
                L = self.total_length, # length of bed [m]
            )

    def _cost(self):
        design_results = self.design_results
        weight = design_results['Weight']
        baseline_purchase_costs = self.baseline_purchase_costs
        if weight == 0:
            baseline_purchase_costs.clear()
            self.parallel.clear()
            return
        baseline_purchase_costs.update(
            self._vessel_purchase_cost(
                weight, design_results['Diameter'], design_results['Length']
            )
        )
        self.parallel['self'] = self.N_columns
        if self.regeneration_fluid:
            baseline_purchase_costs[self.adsorbent] = (
                self.N_columns * self.area * self.column_length
                * self.adsorbent_cost[self.adsorbent]
            )
        else:
            self.ins[2].imass[self.adsorbent] = self.outs[2].imass[self.adsorbent] = self.area * self.column_length / self.cycle_time * self.rho_adsorbent
            self.ins[2].price = self.adsorbent_cost[self.adsorbent] / self.rho_adsorbent

    @staticmethod
    def fit_Freundlich_isotherm(C, q):
        y = np.log(q)
        x = np.log(C)
        m, b = np.polyfit(x, y, 1)
        n = 1 / m
        K = np.exp(b)
        y_ = m * x + b
        q_ = equilibrium_loading_Freundlich_isotherm(C, K, n)
        R2 = np.corrcoef(y, y_)**2
        R2q = np.corrcoef(q, q_)**2
        return {
            'K': K,
            'n': n,
            'parameters': (K, n),
            'linear coefficients': (m, b),
            'linearized data': (x, y),
            'linearized prediction': (x, y_),
            'R2': R2[0, 1],
            'R2-q': R2q[0, 1],
        }
    
    @staticmethod
    def fit_Langmuir_isotherm(C, q):
        y = C / q
        x = C
        m, b = np.polyfit(x, y, 1)
        qmax = 1 / m
        KL = 1 / (b * qmax)
        y_ = m * x + b
        q_ = equilibrium_loading_Langmuir_isotherm(C, KL, qmax)
        R2q = np.corrcoef(q, q_)**2
        R2 = np.corrcoef(y, y_)**2
        return {
            'K': KL,
            'qmax': qmax,
            'parameters': (KL, qmax),
            'linear coefficients': (m, b),
            'linearized data': (x, y),
            'linearized prediction': (x, y_),
            'R2': R2[0, 1],
            'R2-q': R2q[0, 1],
        }
    
    @classmethod
    def plot_isotherm_fit(cls, C, q, method):
        if method == 'Langmuir':
            dct = cls.fit_Langmuir_isotherm(C, q)
        elif method == 'Freundlich':
            dct = cls.fit_Freundlich_isotherm(C, q)
        else:
            raise ValueError("method must be either 'Langmuir' or 'Freundlich'")
        plt.scatter(*dct['linearized data'])
        R2 = dct['R2']
        plt.plot(*dct['linearized prediction'],
                 label=f'R$^2$ = {np.round(R2, 2)}')
        if method == 'Langmuir':
            plt.xlabel(f'C')
            plt.ylabel(f'C / q')
        else:
            plt.xlabel(f'log(C)')
            plt.ylabel(f'log(q)')
    
    @staticmethod
    def fit_solid_mass_transfer_coefficient(t, C, volume, adsorbent, model, args):
        C0 = C[0]
        q = volume * (C0 - C) / adsorbent
        qe = model(C, *args)
        y = np.log((qe - q) / qe)
        x = t
        m = (x * y).sum() / (x * x).sum()
        y_ = m * x
        R2 = np.corrcoef(y, y_)**2
        return {
            'k': -m,
            'linearized data': (x, y),
            'linearized prediction': (x, y_),
            'R2': R2[0, 1],
        }
    
    @staticmethod
    def plot_isotherm_and_mass_transfer_coefficient_fit(
            Ce, qe, t, C, volume, adsorbent, method=None, 
            C_units=None, q_units=None, t_units=None,
        ):
        if C_units is None: C_units = ''
        else: C_units = f' [{format_units(C_units)}]'
        if q_units is None: q_units = ''
        else: q_units = f' [{format_units(q_units)}]'
        if t_units is None: t_units = ''
        else: t_units = f' [{format_units(t_units)}]'
        def plot_fit_isotherm(dct, method):
            plt.scatter(*dct['linearized data'])
            if method == 'Langmuir':
                plt.plot(
                    *dct['linearized prediction'],
                    label=f"R$^2$ = {dct['R2']:.3g}, K={dct['K']:.3g}, q$_{{max}}$={dct['qmax']:.3g}"
                )
                xlabel = 'C {}'
                ylabel = 'C / q {}'
            else:
                plt.plot(
                    *dct['linearized prediction'],
                    label=f"R$^2$ = {dct['R2']:.3g}, K={dct['K']:.3g}, n={dct['n']:.3g}"
                )
                xlabel = 'log(C{})'.format(C_units)
                ylabel = 'log(q{})'.format(q_units)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.legend()
            return dct
        dct_L = bst.AdsorptionColumn.fit_Langmuir_isotherm(Ce, qe)
        dct_F = bst.AdsorptionColumn.fit_Freundlich_isotherm(Ce, qe)
        if method == 'Langmuir' or (method is None and dct_L['R2-q'] > dct_F['R2-q']):
            isotherm_model = 'Langmuir'
            dct = dct_L
        else:
            isotherm_model = 'Freundlich'
            dct = dct_F
        fig, axes = plt.subplots(2, 1)
        ax1, ax2 = axes
        plt.sca(ax1)
        plot_fit_isotherm(dct, isotherm_model)
        model = bst.AdsorptionColumn.isotherm_models[isotherm_model]
        args = dct['parameters']
        dct_k = bst.AdsorptionColumn.fit_solid_mass_transfer_coefficient(
            t, C, volume, adsorbent, model, args
        )
        plt.sca(ax2)
        plt.scatter(*dct_k['linearized data'])
        R2 = dct_k['R2']
        plt.plot(
            *dct_k['linearized prediction'],
            label=f"R$^2$ = {R2:.2g}, k={dct_k['k']:.3g}"
        )
        plt.ylabel('log((qe - q) / q)')
        plt.xlabel(f't{t_units}')
        plt.legend()
        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        return fig, axes, {
            'isotherm_model': isotherm_model, 
            'isotherm_args': args,
            'k': dct_k['k'],
        }

AdsorptionColumn = SingleComponentAdsorptionColumn
