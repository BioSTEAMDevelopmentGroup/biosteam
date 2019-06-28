# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from . import compounds as comps
import numpy as np

Chemical = comps.Chemical
Compound = comps.Compound
inf = np.inf
__all__ = ('Species',)

# TODO: Fix how Dortmund groups are selected in thermo.Chemical. Glycerol defaults to wrong value.

# %% Managing ChEDL thermodynamic properties
setattr_ = object.__setattr__

class Species:
    """Create Species object that contains Compound objects as attributes.

    **Parameters**

        ***IDs:** [str] Strings should be one of the following [-]:
                  * Name, in IUPAC form or common form or a synonym registered in PubChem
                  * InChI name, prefixed by 'InChI=1S/' or 'InChI=1/'
                  * InChI key, prefixed by 'InChIKey='
                  * PubChem CID, prefixed by 'PubChem='
                  * SMILES (prefix with 'SMILES=' to ensure smiles parsing)
                  * CAS number
         
        **Optional**
        
        **cls:** [type Compound] The type of objects that will be created and stored as attributes. Defaults to biosteam.Chemical
        
    """
    units = Compound.units
    _immutable = False
    _compounds = _IDs = _CAS = _Nspecies = None
    
    @classmethod
    def tospecies(cls, compounds, copy=True):
        """Return a Species object from an iterable of Compound objects."""
        if not copy and isinstance(compounds, Species): return compounds
        self = cls.__new__(cls)
        for c in compounds: 
            if isinstance(c, Compound): setattr_(self, c.ID, c)
            else: raise ValueError('compounds must be an iterable of Compound objects')
        return self
    
    def __init__(self, *IDs, cls=None):
        # Sphinx docs look ugly if a class is in the default argument
        if not IDs: return
        if cls is None: cls = Chemical
        
        # Set Compound object attributes
        for n in IDs:
            try:
                compound = cls(n)
            except:
                raise Exception(f"Compound object, '{n}', not defined in data bank")
            setattr_(self, n, compound)

    def _read_only(self):
        if self._immutable: return
        tup = tuple
        compounds = tup(self)
        IDs = tup(i.ID for i in compounds)
        CAS = tup(i.CAS for i in compounds)
        Nspecies = len(IDs)
        index = tup(range(Nspecies))
        me = self.__dict__
        me['_immutable'] = True
        me['_num_compounds'] = tup(zip(index, compounds))
        me['_num_IDs'] = tup(zip(index, IDs))
        me['_compounds'] = compounds
        me['_IDs'] = IDs
        me['_CAS'] = CAS
        me['_Nspecies'] = Nspecies
        me['_MW'] = np.array([i.MW for i in compounds])
        
    def __setattr__(self, ID, compound):
        if self._immutable:
            raise ValueError('cannot alter species object attached to Stream objects')
        elif isinstance(compound, Compound):
            compound.ID = ID
            setattr_(self, ID, compound)
        elif ID == 'IDs':
            setattr_(self, '__dict__', {i: getattr(self, i) for i in compound})
        elif ID == 'CAS':
            index = compound.index if hasattr(compound, 'index') else tuple(compound).index
            setattr_(self, '__dict__', {i.ID: i for i in sorted(self, key=lambda x: index(x.CAS))})
        else:
            raise TypeError('can only set Compound objects as attributes')
    
    def __delattr__(self, ID):
        if self._immutable:
            raise ValueError('cannot alter species object attached to Stream objects')
        else:
            super().__delattr__(ID)
    
    def __len__(self):
        return self._Nspecies or len(self.__dict__)
    
    def __iter__(self):
        return iter(self._compounds or self.__dict__.values())
    
    @property
    def CAS(self):
        return self._CAS or tuple(s.CAS for s in self)
    
    @property
    def IDs(self):
        return self._IDs or tuple(self.__dict__.keys())

    def getprops(self, species_IDs, prop_ID, T, P, phase):
        """Return list of the desired property, prop_ID, for each compound in species_ID."""
        props = []; gat = getattr
        for sID in species_IDs:
            s = gat(self, sID)
            s.T = T
            s.P = P
            s.phase = phase
            props.append(gat(s, prop_ID))
        return props

    def setprops(self, species_IDs, prop_ID, new_values):
        """Set new values to specie property, prop_ID, for given species_IDs."""
        getattr_ = getattr
        for ID in species_IDs: 
            setattr_(getattr_(self, ID), prop_ID, new_values)

    ### Material Properties ###

    # General methods for getting ideal mixture properties
    def _props(self, ID, mol, T, P, phase):
        """Return component property list."""
        out = np.zeros(self._Nspecies)
        getattr_ = getattr
        ic = self._num_compounds.__iter__().__next__
        for m in mol:
            if m: 
                i, c = ic()
                c.P = P
                c.T = T
                c.phase = phase
                prop = getattr_(c, ID)
                out[i] = prop if prop else 0
            else: ic()
        return out

    def _propflows(self, ID, mol, T, P, phase):
        """Return array of component properties * kmol/hr."""
        return self._props(ID, mol, T, P, phase) * mol

    def _propflow(self, ID, mol, T, P, phase):
        """Return flow of component properties."""
        return (self._props(ID, mol, T, P, phase)*mol).sum()

    def _prop(self, ID, mol, T, P, phase):
        """Return molar weighted average property."""
        return (self._props(ID, mol, T, P, phase)*mol).sum().sum()

    def _mixedpropflows(self, ID, molarray, T, P):
        """Return array of component properties * kmol/hr."""
        pf = self._propflows
        return (pf(ID, molarray[0], T, P, 's')
              + pf(ID, molarray[1], T, P, 'l')
              + pf(ID, molarray[2], T, P, 'l')
              + pf(ID, molarray[3], T, P, 'g'))
    
    def _mixedpropflow(self, ID, molarray, T, P):
        """Return flow of component properties."""
        return self._mixedpropflows(ID, molarray, T, P).sum()
    
    def _mixedprop(self, ID, molarray, T, P):
        """Return molar weighted average property."""
        return self._mixedpropflows(ID, molarray, T, P).sum()/molarray.sum()  
    
    ### Equilibrium ###
    
    def _equilibrium_species(self, mol):
        """Return species and indexes of species in equilibrium."""
        species = []; index = []; ic = self._num_compounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.UNIFAC_Dortmund_groups and c.Tb not in (inf, None):
                    species.append(c); index.append(i)
            else: ic()
        return species, index

    def _heavy_species(self, mol):
        """Return species and indexes of heavy species not in equilibrium."""
        species = []; index = []; ic = self._num_compounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.Tb in (inf, None) or not c.UNIFAC_Dortmund_groups:
                    species.append(c); index.append(i)
            else: ic()
        return species, index

    def _light_species(self, mol):
        """Return species and indexes of light species not in equilibrium."""
        species = []; index = []; ic = self._num_compounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.Tb is -inf:
                    species.append(c); index.append(i)
            else: ic()
        return species, index
    
    # Other
    def _reorder(self, values, species_IDs):
        """Return a reordered array with zeros in place of missing species IDs."""
        IDs = self._IDs
        index = IDs.index
        array = np.zeros_like(IDs, dtype=float)
        array[[index(i) for i in species_IDs]] = values
        return array
    
    # Representation
    def __repr__(self):
        return f'<Species: {", ".join(self.IDs)}>'

    
