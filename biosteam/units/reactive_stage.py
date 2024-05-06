# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np

class DelayedFlow:
    __slots__ = ('flows', 'rxn', 'baseline', 'cache')
    
    def __init__(self, flows, rxn, baseline):
        self.flows = flows
        self.rxn = rxn
        self.baseline = baseline
    
    def __call__(self, index):
        try:
            conversion = self.cache
        except:
            self.cache = conversion = self.rxn._conversion(sum(self.flows, self.baseline))
        return self.baseline[index] + conversion[index]
        

class ReactiveStage(bst.Unit):
    _N_outs = _N_ins = 1
    _ins_size_is_fixed = False
    
    def _init(self, T, P, phase, rxn):
        self.T = T
        self.P = P
        self.phase = phase
        self.rxn = rxn
        
    def _run(self):
        outlet, = self.outs
        outlet.T = self.T
        outlet.P = self.P
        outlet.mol = sum([i.mol for i in self.ins])
        outlet.phase = self.phase
        self.rxn(outlet)
    
    def _create_material_balance_equations(self):
        inlets = self.ins
        product, = self.outs
        equations = []
        n = self.chemicals.size
        ones = np.ones(n)
        minus_ones = -ones
        fresh_inlets = [i for i in inlets if i.isfeed() and not i.equations]
        process_inlets = [i for i in inlets if not i.isfeed() or i.equations]
        
        nonlimiting_players = []
        rxn = self.rxn
        if isinstance(rxn, (bst.Rxn, bst.RxnI)):
            index = rxn._X_index[1] if rxn.phases else rxn._X_index
            for i in rxn._stoichiometry.nonzero_keys():
                if i == index: continue
                nonlimiting_players.append(i)
            nonlimiting_set = set(nonlimiting_players)
        else:
            # TODO: implement case with reaction system / parallel reaction / series reaction
            # TODO: implement case with chemicals with multiple roles as limiting reactants/products/co-reactants
            raise NotImplementedError('only single reaction objects work (for now)')
        # Overall flows
        eq_overall = {}
        flows = [i.imol.data for i in process_inlets]
        predetermined_flow = sum([i.mol for i in fresh_inlets])
        eq_overall[product] = ones
        for i in process_inlets:
            inlet_coeffs = minus_ones.copy()
            inlet_coeffs[index] += rxn.X
            eq_overall[i] = inlet_coeffs
        rhs = predetermined_flow.to_array(dtype=object)
        delayed_flow = DelayedFlow(flows, rxn, predetermined_flow)
        for i in nonlimiting_set: rhs[i] = delayed_flow
        rhs[index] *= (1 - rxn.X)
        equations.append(
            (eq_overall, rhs)
        )
        return equations
    
    def _create_linear_equations(self, variable):
        # list[dict[Unit|Stream, float]]
        if variable == 'material':
            return self._create_material_balance_equations()
        else:
            eqs = []
        return eqs
    
    def _update_decoupled_variable(self, variable, value):
        pass
    
if __name__ == '__main__':
    bst.settings.set_thermo(['N2', 'H2', 'NH3'])
    bst.settings.mixture.include_excess_energies = True
    feed = bst.Stream(N2=1, H2=3, total_flow=1000)
    recycle = bst.Stream()
    rxstage = ReactiveStage(
        ins=[feed, recycle], 
        rxn=bst.Rxn('N2 + H2 -> NH3', reactant='N2', X=0.2, correct_atomic_balance=True),
        T=425, P=275 * 101325, phase='g',
    )
    # flash = bst.StageEquilibrium(
    #     ins=[rxstage-0], outs=[recycle, 'ammonia'], phases=('g', 'l'),
    #     T=30 + 273.15, P=250 * 101325
    # )
    flash = bst.Separator(
        ins=[rxstage-0], outs=[recycle, 'ammonia'], phases=('g', 'l'),
        T=30 + 273.15, P=250 * 101325, split=dict(NH3=0, N2=0.999, H2=0.999),
    )
    sm = bst.System.from_units(units=[rxstage, flash], algorithm='Sequential modular')
    sm.simulate()
    sm.show()
    recycle.empty()
    po = bst.System.from_units(units=[rxstage, flash], algorithm='Phenomena oriented')
    po.simulate()
    po.show()
    