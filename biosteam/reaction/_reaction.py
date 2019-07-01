# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:23:52 2019

@author: yoelr
"""
from . import _parse as prs
from .. import Stream
from .._exceptions import UndefinedCompound
import numpy as np

__all__ = ('Reaction', 'ParallelReaction', 'SeriesReaction', 'Rxn', 'SRxn', 'PRxn')

def stoi2str(stoi, species):
    return f"{prs.arr2str(stoi, species)}"

class Reaction:
    """Create a Reaction object which defines a stoichiometric reaction and conversion.
    
    **Parameters**
    
        **reaction:** [str] Stoichiometric equation and conversion of a reactant written as:
            i1 R1 + ... + in Rn -> j1 P1 + ... + jm Pm

        **reactant:** [str] ID of reactant compound.
        
        **X:** [float] Reactant conversion (fraction).
        
        **species:** [Species] Defaults to Stream.species.
        
    **Examples**
    
        >>> import biosteam as bst    
        >>> import biosteam.reaction as rn
        >>> bst.Stream.species = sp = bst.Species('H2O', 'H2', 'O2')
        >>> srn = rn.Reaction('2H2O -> 2H2 + O2', reactant='H2O', X=0.7)
        >>> srn
        Reaction('H2O -> H2 + 0.5 O2', reactant='H2O', X=0.7)
        >>> feed = bst.Stream('feed', H2O=200)
        >>> feed.mol[:] += srn(feed.mol) # Call to run reaction on molar flow
        >>> feed
        Stream: feed
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  60
                         H2   140
                         O2   70
        
        Notice how 70% of water was converted to product.
    
    """
    __slots__ = ('_species', '_Xindex', '_stoi', 'X')
    def __init__(self, reaction, reactant, X, species=None):
        if not species: species = Stream.species
        self._species = species
        self._stoi = prs.str2arr(reaction, species)
        self.reactant = reactant
        self.X = X #: [float] Reactant conversion
    
    def __call__(self, material):
        return material[self._Xindex]*self.X*self._stoi
    
    @property
    def species(self):
        """Species corresponing to each entry in the stoichiometry array."""
        return self._species
    @property
    def stoichiometry(self):
        """[array] Stoichiometry coefficients."""
        return self._stoi
    
    @property
    def reactant(self):
        """[str] Reactant associated to conversion."""
        return self._species._IDs[self._Xindex]
    @reactant.setter
    def reactant(self, ID):
        try: self._Xindex = self._species._indexdct[ID]
        except KeyError: raise UndefinedCompound(ID)
        self._stoi *= 1/-(self._stoi[self._Xindex])
    
    def __repr__(self):
        return f"{type(self).__name__}('{stoi2str(self._stoi, self._species)}', reactant='{self.reactant}', X={self.X:.3g})"
    
    def show(self):
        stoichiometry = stoi2str(self._stoi, self._species)
        print(f"{type(self).__name__}:\n"
             +f" {stoichiometry}\n"
             +f" {self.X:.2%} {self.reactant} conversion")
    _ipython_display_ = show

Rxn = Reaction

class ReactionSet:
    """Abstract class for a set of reactions."""
    __slots__ = ('_stoi', 'X', '_Xindex', '_species')
    def __init__(self, reactions):
        if not reactions:
            raise ValueError('reactions must not be empty')
        species = {i._species for i in reactions}
        if len(species) != 1:
            raise ValueError('all reactions must have the same species')
        self._stoi = np.array([i._stoi for i in reactions])
        self.X = np.array([i.X for i in reactions])
        self._Xindex = np.array([i._Xindex for i in reactions])
        self._species = reactions[0]._species
    
    @property
    def species(self):
        """Species corresponing to each entry in the stoichiometry array."""
        return self._species
    @property
    def stoichiometry(self):
        """[array] Stoichiometry coefficients."""
        return self._stoi
    
    @property
    def reactant(self):
        """[str] Reactant associated to conversion."""
        IDs = self._species._IDs
        return tuple(IDs[i] for i in self._Xindex)
    
    def __repr__(self):
        return f"<{type(self).__name__}: {', '.join(set(self.reactant))}>"
    
    def show(self):
        outs = f"{type(self).__name__}:"
        species = self._species
        rxns = [stoi2str(i, species) for i in self._stoi]
        maxrxnlen = max([13, *[len(i) for i in rxns]]) + 2
        cmps = self.reactant
        maxcmplen = max([8, *[len(i) for i in cmps]]) + 2
        Xs = self.X
        outs += "\n stoichiometry" + " "*(maxrxnlen - 13) + "reactant" + " "*(maxcmplen - 8) + '  X[%]'
        for rxn, cmp, X in zip(rxns, cmps, Xs):
            rxn_spaces = " "*(maxrxnlen - len(rxn))
            cmp_spaces = " "*(maxcmplen - len(cmp))
            outs += f"\n {rxn}{rxn_spaces}{cmp}{cmp_spaces}{X*100: >6.2f}"
        print(outs)
    _ipython_display_ = show
        
class ParallelReaction(ReactionSet):
    __slots__ = ()
    
    def __call__(self, material):
        return (material[self._Xindex]*self.X*self._stoi).sum()

    @property
    def X_net(self):
        X_net = {}
        for i, j in zip(self.reactant, self.X):
            if i in X_net:
                X_net[i] += j
            else:
                X_net[i] = j
        return prs.dct2arr(X_net, self.species)
    

class SeriesReaction(ReactionSet):
    __slots__ = ()
    
    def __call__(self, material):
        for i, j, k in (self._Xindex, self.X, self._stoi):
            material = material + material[self._Xindex]*self.X*self._stoi
        return material

    @property
    def X_net(self):
        X_net = {}
        for i, j in zip(self.reactant, self.X):
            if i in X_net:
                X_net[i] *= j
            else:
                X_net[i] = j
        return prs.dct2arr(X_net, self.species)

PRxn = ParallelReaction
SRxn = SeriesReaction




    