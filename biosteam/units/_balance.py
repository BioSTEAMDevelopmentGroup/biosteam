# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 17:35:28 2018

@author: yoelr
"""
import numpy as np
from ..utils.piping import StreamSequence
from .. import Unit, PowerUtility

__all__ = ('MassBalance',)

# %% Mass Balance Unit

class MassBalance(Unit):
    """Create a Unit object that changes net input flow rates to satisfy output flow rates. This calculation is based on mass balance equations for specified IDs. 

    Parameters
    ----------
    ins :
        Inlet streams.
    outs : 
        Outlet streams
    chemical_IDs : tuple[str]
        Chemicals that will be used to solve mass balance linear equations.
        The number of chemicals must be same as the number of input streams varied.
    streams : tuple[int]
        Indices of input streams that can vary in net flow rate
    is_exact=True : bool, optional
        True if exact flow rate solution is required for the specified IDs.
    balance='flow' : {'flow', 'composition'}, optional
          * 'flow': Satisfy output flow rates
          * 'composition': Satisfy net output molar composition

    Examples
    --------
    MassBalance are Unit objects that serve to alter flow rates of selected
    chemicals and input streams to satisfy the mass balance.
    The example below uses the MassBalance object to satisfy the target
    flow rate feeding the mixer M1:
    
    >>> from biosteam import System
    >>> from biosteam.units import Mixer, Splitter, StorageTank, Pump, Flash, MassBalance
    >>> from thermosteam import Chemicals, Stream, settings
    >>> chemicals = Chemicals(['Water', 'Ethanol'])
    >>> settings.set_thermo(chemicals)
    >>> water = Stream('water',
    ...                Water=40,
    ...                units='lb/s',
    ...                T=350, P=101325)
    >>> ethanol = Stream('ethanol',
    ...                  Ethanol=190, Water=30,
    ...                  T=300, P=101325)
    >>> target = Stream('target', flow=[50, 50])
    >>> T1 = StorageTank('T1')
    >>> T2 = StorageTank('T2')
    >>> P1 = Pump('P1', P=101325)
    >>> P2 = Pump('P2', P=101325)
    >>> M1 = Mixer('M1', outs='s1')
    >>> S1 = Splitter('S1', outs=('s2', 's3'), split=0.5)
    >>> F1 = Flash('F1', outs=('s4', 's5'), V=0.5, P =101325)
    >>> # Connect units
    >>> water-T1-P1
    <Pump: P1>
    >>> ethanol-T2-P2
    <Pump: P2>
    >>> [P1-0, P2-0, S1-0]-M1-F1-1-S1
    <Splitter: S1>
    >>> MB1 = MassBalance('MB1', streams=[0,1],
    ...                   chemical_IDs=['Ethanol', 'Water'],
    ...                   outs=target,
    ...                   ins=(water, ethanol, S1-0))
    >>> mixSys = System('mixSys',
    ...                  recycle=S1-0,
    ...                  network=(MB1, T1, P1, T2, P2, M1, F1, S1))
    >>> # Make diagram to view system
    >>> # mixSys.diagram()
    >>> mixSys.simulate();
    >>> MB1.show()
    MassBalance: MB1
    ins...
    [0] water
        phase: 'l', T: 350 K, P: 101325 Pa
        flow (kmol/hr): Water  28.3
    [1] ethanol
        phase: 'l', T: 300 K, P: 101325 Pa
        flow (kmol/hr): Water    6.37
                        Ethanol  40.3
    [2] s2  from  Splitter-S1
        phase: 'l', T: 353.88 K, P: 101325 Pa
        flow (kmol/hr): Water    15.3
                        Ethanol  9.65
    outs...
    [0] target
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water    50
                        Ethanol  50
    
    """
    line = 'Balance'
    results = None
    _has_cost = False
    _N_ins = 2
    heat_utilities = ()
    power_utility = None

    def __init__(self, ID='', outs=(), ins=None,
                 chemical_IDs=None, streams=None,
                 is_exact=True, balance='flow', thermo=None):
        self.ID = ID
        thermo = self._load_thermo(thermo)
        self._ins = StreamSequence(self._N_ins, ins, thermo, fixed_size=False)
        self._outs = StreamSequence(self._N_outs, outs, thermo, fixed_size=False)
        self._power_utility = PowerUtility()
        self.streams = streams
        self.chemical_IDs = chemical_IDs
        self.is_exact = is_exact
        self.balance = balance
        
    def _run(self):
        """Solve mass balance by iteration."""
        # SOLVING BY ITERATION TAKES 15 LOOPS FOR 2 STREAMS
        # SOLVING BY LEAST-SQUARES TAKES 40 LOOPS
        stream_index = self.streams
        balance = self.balance
        solver = np.linalg.solve if self.is_exact else np.linalg.lstsq
        ins = self.ins

        # Set up constant and variable streams
        vary = [ins[i] for i in stream_index]  # Streams to vary in mass balance
        constant = [i for i in ins if (i and i not in vary)]  # Constant streams
        index = self.chemicals.get_index(self.chemical_IDs)

        if balance == 'flow':
            # Perform the following calculation: Ax = b = f - g
            # Where:
            #    A = flow rate array
            #    x = factors
            #    b = target flow rates
            #    f = output flow rates
            #    g = constant input flow rates

            # Solve linear equations for mass balance
            A = np.array([s.mol for s in vary]).transpose()[index, :]
            f = self.mol_out[index]
            g = sum([s.mol[index] for s in constant])
            b = f - g
            x = solver(A, b)

            # Set flow rates for input streams
            for factor, s in zip(x, vary):
                s.mol[:] = s.mol * factor

        elif balance == 'fraction':
            # Perform the following calculation:
            # Ax = b
            #    = sum( A_ * x_guess + g_ )f - g
            #    = A_ * x_guess * f - O
            # O  = sum(g_)*f - g
            # Where:
            # A_ is flow array for all species
            # g_ is constant flows for all species
            # Same variable definitions as in 'flow'

            # Set all variables
            A_ = np.array([s.mol for s in vary]).transpose()
            A = np.array([s.mol for s in vary]).transpose()[index, :]
            f = self.z_mol_out[index]
            g_ = sum([s.mol for s in constant])
            g = g_[index]
            O = sum(g_) * f - g

            # Solve by iteration
            x_guess = np.ones_like(index)
            not_converged = True
            while not_converged:
                # Solve linear equations for mass balance
                b = (A_ * x_guess).sum()*f + O
                x_new = solver(A, b)
                not_converged = sum(((x_new - x_guess)/x_new)**2) > 0.0001
                x_guess = x_new

            # Set flow rates for input streams
            for factor, s in zip(x_new, vary):
                s.mol = s.mol * factor
        
        else:
            raise ValueError( "balance type must be one of the following: 'flow', 'composition'")

# Balance
graphics = MassBalance._graphics
graphics.node['shape'] = 'note'
graphics.node['fillcolor'] = '#F0F0F0'
graphics.in_system = False
del graphics


# %% Energy Balance Unit

# class EnergyBalance(Unit):
#     """Create a Unit object that changes a stream's temperature, flow rate, or vapor fraction to satisfy energy balance.

#     **Parameters**

#         **index:** [int] Index of stream that can vary in temperature, flow rate, or vapor fraction.
        
#         **Type:** [str] Should be one of the following
#             * 'T': Vary temperature of output stream
#             * 'F': Vary flow rate of input/output stream
#             * 'V': Vary vapor fraction of output stream
        
#         **Qin:** *[float]* Additional energy input.
        
#     .. Note:: This is not a mixer, input streams and output streams should match flow rates.

#     """
#     _kwargs = {'index': None,
#                'Type': 'T',
#                'Qin': 0}
#     line = 'Balance'
#     _has_cost = False
#     _graphics = MassBalance._graphics
#     _init_ins = MassBalance._init_ins
#     _init_outs = MassBalance._init_outs
    
#     def _run(self):        # Get arguments
#         ins = self.ins.copy()
#         outs = self.outs.copy()
#         kwargs = self._kwargs
#         index = kwargs['index']
#         Type = kwargs['Type']
#         Qin = kwargs['Qin']
        
#         # Pop out required streams
#         if Type == 'F':
#             s_in = ins.pop(index)
#             s_out = outs.pop(index)
#         else:
#             s = outs.pop(index)
        
#         # Find required enthalpy
#         H_in = sum(i.H for i in ins) + Qin
#         H_out = sum(o.H for o in outs)
#         H_s = H_out - H_in
        
#         # Set enthalpy
#         if Type == 'T':
#             s.H = -H_s
#         elif Type == 'V':
#             s.enable_phases()
#             s.VLE(Qin=s.H - H_s)
#         elif Type == 'F':
#             s.mol *= (s_out.H - s_in.H)/H_s
#         else:
#             raise ValueError(f"Type must be 'T', 'V' or 'F', not '{Type}'")
            
        
        