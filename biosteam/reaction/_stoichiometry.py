# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:23:52 2019

@author: yoelr
"""
from . import _parse as prs
from .. import Stream
from .._exceptions import UndefinedCompound

__all__ = ('Stoichiometry',)

def stoich2str(sr):
    return f"{prs.arr2str(sr._coeff, sr._species)}"

class Stoichiometry:
    """Create a Stoichiometry object which defines a stoichiometric reaction and conversion.
    
    **Parameters**
    
        **reaction:** [str] Stoichiometric equation and conversion of a reactant written as:
            i1 R1 + ... + in Rn -> j1 P1 + ... + jm Pm; X Ri

        **species:** [Species] Defaults to Stream.species.
        
    **Examples**
    
        >>> import biosteam as bst    
        >>> import biosteam.reaction as rn
        >>> bst.Stream.species = sp = bst.Species('H2O', 'H2', 'O2')
        >>> srn = rn.Stoichiometry('2H2O -> 2H2 + O2', reactant='H2O', X=0.7)
        >>> srn
        Stoichiometry('H2O -> H2 + 0.5 O2', reactant='H2O', X=0.7)
        >>> feed = bst.Stream('feed', H2O=200)
        >>> srn(feed.mol) # Call to run reaction on molar flow
        >>> feed
        Stream: feed
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  60
                         H2   140
                         O2   70
        
        Notice how 70% of water was converted to product.
    
    """
    __slots__ = ('_species', '_Xindex', '_coeff', 'X')
    def __init__(self, reaction, reactant, X, species=None):
        reaction = reaction.replace(' ', '')
        if not species: species = Stream.species
        self._species = species
        self._coeff = prs.str2arr(reaction, species)
        self.reactant = reactant
        self.X = X
    
    def __call__(self, mol):
        mol += mol[self._Xindex]*self.X*self._coeff
    
    @property
    def species(self):
        return self._species
    @property
    def coeff(self):
        return self._coeff
    
    @property
    def reactant(self):
        """[str] Reactant associated to conversion."""
        return self._species._IDs[self._Xindex]
    @reactant.setter
    def reactant(self, ID):
        try: self._Xindex = self._species._indexdct[ID]
        except KeyError: raise UndefinedCompound(ID)
        self._coeff *= 1/-(self._coeff[self._Xindex])
    
    def __repr__(self):
        return f"{type(self).__name__}('{stoich2str(self)}', reactant='{self.reactant}', X={self.X:.3g})"

# class DynamicFlow:
#     """Create a DynamicFlow object which defines a flow rate requirement.
    
#     **Parameters**
        
#         ****compound_ratios:** Pairs of compound ID and ratio of compound to reactant.

#         **IDs:** tuple[str] Compound IDs. Defaults to Stream.species.IDs
        
#     **Examples**
    
#         >>> from biosteam import *
#         >>> from lipidcane.species import biodiesel_species
#         >>> Stream.species = biodiesel_species
#         >>> feed = Stream('Catalyst')
#         >>> dynflow = reaction.DynamicFlow(Methanol=6, NaOCH3=0.06384)
#         >>> dynflow
#         DynamicFlow(Methanol=6, NaOCH3=0.06384)
#         >>> dynflow(100) # Call to return flow
#         >>> np.array([])
    
    
#     """
#     __slots__ = ('_reactant', '_ratios', '_IDs')
#     def __init__(self, IDs=None, **compound_ratios):
#         self._IDs = _IDs = _get_IDs(IDs)
#         self._ratios = rn_dct2arr(compound_ratios, _IDs)
    
#     def __call__(self, flow):
#         return flow * self._ratios
    
#     @property
#     def ratios(self):
#         return self._ratios
#     @property
#     def IDs(self):
#         return self._IDs
    
    