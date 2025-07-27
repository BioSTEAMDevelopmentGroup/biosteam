# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module contains unit operations for splitting flows.

.. contents:: :local:
    
.. autoclass:: biosteam.units.splitting.Splitter
.. autoclass:: biosteam.units.splitting.PhaseSplitter 
.. autoclass:: biosteam.units.splitting.MockSplitter
.. autoclass:: biosteam.units.splitting.ReversedSplitter

"""
from .. import Unit
from thermosteam._graphics import splitter_graphics
from thermosteam import separations
import biosteam as bst
import numpy as np
from thermosteam import VariableNode

__all__ = ('Splitter', 'PhaseSplitter', 'FakeSplitter', 'MockSplitter',
           'ReversedSplitter', 'Separator')

class Splitter(Unit):
    """
    Create a splitter that separates mixed streams based on splits.

    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Split stream
        * [1] Remainder stream    
    split : Should be one of the following
        * [float] The fraction of net feed in the 0th outlet stream
        * [array_like] Componentwise split of feed to 0th outlet stream
        * [dict] ID-split pairs of feed to 0th outlet stream
    order=None : Iterable[str], defaults to biosteam.settings.chemicals.IDs
        Chemical order of split.
    
    Examples
    --------
    Create a Splitter object with an ID, a feed stream, two outlet streams,
    and an overall split:
        
    >>> from biosteam import units, settings, Stream
    >>> settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = Stream('feed', Water=20, Ethanol=10, T=340)
    >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'), split=0.1)
    >>> S1.simulate()
    >>> S1.show()
    Splitter: S1
    ins...
    [0] feed
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    outs...
    [0] top
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    2
                        Ethanol  1
    [1] bot
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    18
                        Ethanol  9
      
    Create a Splitter object, but this time with a componentwise split
    using a dictionary:
    
    >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
    ...                     split={'Water': 0.1, 'Ethanol': 0.99})
    >>> S1.simulate()
    >>> S1.show()
    Splitter: S1
    ins...
    [0] feed
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    outs...
    [0] top
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    2
                        Ethanol  9.9
    [1] bot
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    18
                        Ethanol  0.1
                           
    Create a Splitter object using componentwise split, but this time specify the order:
    
    >>> S1 = units.Splitter('S1', ins=feed, outs=('top', 'bot'),
    ...                     order=('Ethanol', 'Water'),
    ...                     split=(0.99, 0.10))
    >>> S1.simulate()
    >>> S1.show()
    Splitter: S1
    ins...
    [0] feed
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    20
                        Ethanol  10
    outs...
    [0] top
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    2
                        Ethanol  9.9
    [1] bot
        phase: 'l', T: 340 K, P: 101325 Pa
        flow (kmol/hr): Water    18
                        Ethanol  0.1

    Splits can also be altered after creating the splitter:
        
    >>> S1.split = 0.5
    >>> S1.isplit.show()
    SplitIndexer:
     Water    0.5
     Ethanol  0.5
     
    >>> S1.isplit['Water'] = 1.0
    >>> S1.isplit.show()
    SplitIndexer:
     Water    1
     Ethanol  0.5
     
    >>> S1.split = [0.9, 0.8]
    >>> S1.isplit.show()
    SplitIndexer:
     Water    0.9
     Ethanol  0.8

    """
    _N_outs = 2
    _graphics = splitter_graphics
    
    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[SparseArray] Componentwise split of feed to 0th outlet stream."""
        return self._isplit.data
    @split.setter
    def split(self, values):
        split = self.split
        if split is not values:
            split[:] = values
    
    def _init(self, split, order=None):
        self._isplit = self.thermo.chemicals.isplit(split, order)
        
    def _run(self):
        feed = self._ins[0]
        isplit = self._isplit
        if isplit.chemicals is not feed.chemicals: self._reset_thermo(feed._thermo)
        feed.split_to(*self.outs, isplit.data)
        
    # def _create_material_balance_equations(self, composition_sensitive):
    #     fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
    #     top, bottom = self.outs
    #     ones = np.ones(self.chemicals.size)
    #     minus_ones = -ones
    #     zeros = np.zeros(self.chemicals.size)
        
    #     # Overall flows
    #     eq_overall = {}
    #     for i in self.outs: 
    #         eq_overall[i] = ones
    #     for i in process_inlets:
    #         if i in eq_overall: del eq_overall[i]
    #         else: eq_overall[i] = minus_ones
    #     equations.append(
    #         (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
    #     )
        
    #     # Top and bottom flows
    #     eq_outs = {}
    #     split = self.split
    #     minus_split = -split
    #     for i in process_inlets: eq_outs[i] = minus_split
    #     rhs = split * sum([i.mol for i in fresh_inlets], zeros)
    #     eq_outs[top] = ones
    #     equations.append(
    #         (eq_outs, rhs)
    #     )
    #     return equations
    
    # def _get_energy_departure_coefficient(self, stream):
    #     coeff = -stream.C
    #     return (self, coeff)
    
    # def _create_energy_departure_equations(self):
    #     coeff = {self: self.ins[0].C}
    #     self.ins[0]._update_energy_departure_coefficient(coeff)
    #     return [(coeff, self.H_in - self.H_out)]
    
    # def _update_energy_variable(self, departure):
    #     for i in self.outs: i.T += departure


class PhaseSplitter(Unit):
    """
    Create a PhaseSplitter object that splits the feed to outlets by phase.
    
    Parameters
    ----------
    ins : 
        Feed.
    outs : 
        Outlets.
        
    Notes
    -----
    Phases allocate to outlets in alphabetical order. For example,
    if the feed.phases is 'gls' (i.e. gas, liquid, and solid), the phases
    of the outlets will be 'g', 'l', and 's'.
        
    Examples
    --------
    >>> import biosteam as bst
    >>> bst.settings.set_thermo(['Water', 'Ethanol'], cache=True)
    >>> feed = bst.Stream('feed', Water=10, Ethanol=10)
    >>> feed.vle(V=0.5, P=101325)
    >>> s1 = bst.Stream('s1')
    >>> s2 = bst.Stream('s2')
    >>> PS = bst.PhaseSplitter('PS', feed, [s1, s2])
    >>> PS.simulate()
    >>> PS.show()
    PhaseSplitter: PS
    ins...
    [0] feed
        phases: ('g', 'l'), T: 353.94 K, P: 101325 Pa
        flow (kmol/hr): (g) Water    3.87
                            Ethanol  6.13
                        (l) Water    6.13
                            Ethanol  3.87
    outs...
    [0] s1
        phase: 'g', T: 353.94 K, P: 101325 Pa
        flow (kmol/hr): Water    3.87
                        Ethanol  6.13
    [1] s2
        phase: 'l', T: 353.94 K, P: 101325 Pa
        flow (kmol/hr): Water    6.13
                        Ethanol  3.87
    
    """
    _N_ins = 1
    _N_outs = 2
    _graphics = splitter_graphics
    
    def _run(self):
        separations.phase_split(*self.ins, self.outs)


class MockSplitter(Unit):
    """
    Create a MockSplitter object that does nothing when simulated.
    """
    _graphics = Splitter._graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    
    def _run(self): pass

MockSplitter.line = 'Splitter'
FakeSplitter = MockSplitter    

class ReversedSplitter(Unit):
    """
    Create a splitter that, when simulated, sets the inlet stream based 
    on outlet streams. Must have only one input stream. The outlet streams will
    have the same temperature, pressure and phase as the inlet.
    
    """
    _graphics = Splitter._graphics
    _N_ins = 1
    _N_outs = 2
    _outs_size_is_fixed = False
    power_utility = None
    heat_utilities = ()
    results = None
    
    def _run(self):
        inlet, = self.ins
        outlets = self.outs
        reversed_split(inlet, outlets)


def reversed_split(inlet, outlets):
    inlet.mol[:] = sum([i.mol for i in outlets])
    T = inlet.T
    P = inlet.P
    phase = inlet.phase
    for out in outlets:
        out.T = T
        out.P = P
        out.phase = phase 


class Separator(Unit):
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    
    @property
    def equation_node_names(self): 
        material_balances = (
            'overall_material_balance_node', 
            'separation_material_balance_node',
        )
        if self.T is None:
            return (
                *material_balances,
                'energy_balance_node',
            )
        else:
            return material_balances
    
    def initialize_overall_material_balance_node(self):
        self.overall_material_balance_node.set_equations(
            inputs=[j for i in self.ins if (j:=i.F_node)],
            outputs=[i.F_node for i in self.outs],
        )
    
    def initialize_separation_material_balance_node(self):
        self.separation_material_balance_node.set_equations(
            outputs=[i.F_node for i in self.outs],
            inputs=[],
        )
        
    def initialize_energy_balance_node(self):
        self.energy_balance_node.set_equations(
            inputs=(
                self.T_node, 
                *[i.T_node for i in (*self.ins, *self.outs)],
                *[i.F_node for i in (*self.ins, *self.outs)],
                *[j for i in self.ins if (j:=i.E_node)]
            ),
            outputs=[j for i in self.outs if (j:=i.E_node)],
        )
        
    @property
    def T_node(self):
        if not self.T:
            if hasattr(self, '_T_node'): return self._T_node
            self._T_node = var = VariableNode(f"{self.node_tag}.T", lambda: self.outs[0].T)
            return var 
    
    def get_E_node(self, stream):
        return self.T_node
    
    @property
    def E_node(self):
        return self.T_node
    
    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[SparseArray] Componentwise split of feed to 0th outlet stream."""
        return self._isplit.data
    @split.setter
    def split(self, values):
        split = self.split
        if split is not values:
            split[:] = values
    
    def _init(self, split, order=None, T=None, P=None, phases=None):
        self._isplit = self.thermo.chemicals.isplit(split, order)
        self.T = T
        self.P = P
        self.phases = phases
        
    def _mass_and_energy_balance_specifications(self):
        isplit = self._isplit
        specs = [
            (i.ID + ' split', j * 100, '%') for i, j in zip(isplit.chemicals, isplit.data)
        ]
        if self.T is not None: 
            specs.append(
                ('T', self.T, 'K')
            )
        if self.phases is not None:
            specs.append(
                ('Phases', self.phases, '-')
            )
        if self.P is not None:
            specs.append(
                ('P', self.P, 'Pa')
            )
        return 'Separator', specs
        
    def _run(self):
        ins = self._ins
        top, bottom = self._outs
        top.mix_from(ins)
        top.split_to(*self.outs, self._isplit.data)
        if self.T:
            top.T = bottom.T = self.T
        if self.phases:
            top.phase, bottom.phase = self.phases
        if self.P:
            top.P = bottom.P = self.P
        
    def _create_material_balance_equations(self, composition_sensitive):
        fresh_inlets, process_inlets, equations = self._begin_equations(composition_sensitive)
        top, bottom = self.outs
        ones = np.ones(self.chemicals.size)
        minus_ones = -ones
        zeros = np.zeros(self.chemicals.size)
        
        # Overall flows
        eq_overall = {}
        for i in self.outs: 
            eq_overall[i] = ones
        for i in process_inlets:
            if i in eq_overall: del eq_overall[i]
            else: eq_overall[i] = minus_ones
        equations.append(
            (eq_overall, sum([i.mol for i in fresh_inlets], zeros))
        )
        
        # Top and bottom flows
        eq_outs = {}
        split = self.split
        minus_split = -split
        for i in process_inlets: eq_outs[i] = minus_split
        rhs = split * sum([i.mol for i in fresh_inlets], zeros)
        eq_outs[top] = ones
        equations.append(
            (eq_outs, rhs)
        )
        return equations
    
    def _get_energy_departure_coefficient(self, stream):
        if self.T is not None: return
        coeff = -stream.C
        return (self, coeff)
    
    def _create_energy_departure_equations(self):
        if self.T is not None: return []
        coeff = {self: self.ins[0].C}
        self.ins[0]._update_energy_departure_coefficient(coeff)
        return [(coeff, self.H_in - self.H_out)]
    
    def _update_energy_variable(self, departure):
        for i in self.outs: i.T += departure
    
