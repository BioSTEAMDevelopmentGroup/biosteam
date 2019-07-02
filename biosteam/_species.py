# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from . import compounds as cmps
from ._exceptions import UndefinedCompound
import numpy as np

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
    _immutable = False
    _compounds = _IDs = _Nspecies = None
    
    @classmethod
    def tospecies(cls, compounds, copy=True):
        """Return a Species object from an iterable of Compound objects."""
        if not copy and isinstance(compounds, Species): return compounds
        self = cls.__new__(cls)
        for c in compounds: 
            if isinstance(c, cmps.Compound): setattr_(self, c.ID, c)
            else: raise ValueError('compounds must be an iterable of Compound objects')
        return self
    
    def __dir__(self):
        return [i for i in super().__dir__() if i.isalnum()]
    
    def __init__(self, *IDs, cls=None):
        # Sphinx docs look ugly if a class is in the default argument
        if not IDs: return
        if cls is None: cls = cmps.Chemical
        
        # Set Compound object attributes
        for n in IDs:
            try:
                compound = cls(n)
            except:
                raise LookupError(f"'{n}', not defined in data bank")
            setattr_(self, n, compound)

    def set_synonym(self, ID, synonym):
        compound = getattr(self, ID)
        dct = self.__dict__
        if synonym in dct and dct[synonym] is not compound:
            raise ValueError(f"synonym '{synonym}' already in use by {repr(dct[synonym])}")
        else:
            try:
                self._indexdct[synonym] = self._indexdct[ID]
            except AttributeError:
                raise RuntimeError("'{type(self).__name__}' object must be read-only to set synonyms")
            dct[synonym] = compound
            

    def kwarray(self, **data):
        """Return an array with entries that correspond to compound IDs.
        
        **Parameters**
        
            ****data:** ID-value pair.
            
        **Examples**
        
            >>> from biosteam import *
            >>> Stream.species = Species('Water', 'Ethanol')
            >>> Stream.species.kwarray(Water=2)
            array([2., 0.])
        
        """
        return self.array(data, [*data.values()])
    
    def array(self, IDs, data):
        """Return an array with entries that correspond to species IDs.
        
        **Parameters**
        
            **IDs:** [iterable] compound IDs.
            
            **data:** [array_like] Data corresponding to IDs.
            
        **Examples**
        
            >>> from biosteam import *
            >>> Stream.species = Species('Water', 'Ethanol')
            >>> Stream.species.array(['Water'], [2])
            array([2., 0.])
        
        """
        array = np.zeros(self._Nspecies)
        array[self.indices(IDs)] = data
        return array

    def indices(self, IDs):
        """Return flow indices of specified species.

        **Parameters**

             **IDs:** [iterable] species IDs or CAS numbers.

        Example

        Indices by ID:
        
        >>> from biosteam import *
        >>> Stream.species = Species(['Ethanol', 'Water'])
        >>> Stream.species.indices(['Water', 'Ethanol'])
        [1, 0]

        Indices by CAS number:
        
        >>> Stream.species.indices(['7732-18-5', '64-17-5']):
        [1, 0]

        """
        try:
            dct = self._indexdct
            return [dct[i] for i in IDs]
        except AttributeError:
            raise RuntimeError("'{type(self).__name__}' object must be read-only to get indices")
        except KeyError:
            for i in IDs:
                if i not in dct: raise UndefinedCompound(i)
                

    def read_only(self):
        """Make object read-only (attributes will be immutable)."""
        if self._immutable: return
        tup = tuple
        compounds = tup(self.__dict__.values())
        IDs = tup(i for i in self.__dict__)
        CAS = tup(i.CAS for i in compounds)
        Nspecies = len(IDs)
        index = tup(range(Nspecies))
        me = self.__dict__
        for i in compounds: me[i.CAS] = i
        me['_indexdct'] = dict((*zip(CAS, index), *zip(IDs, index)))
        me['_immutable'] = True
        me['_numcompounds'] = tup(zip(index, compounds))
        me['_numIDs'] = tup(zip(index, IDs))
        me['_compounds'] = compounds
        me['_IDs'] = IDs
        me['_Nspecies'] = Nspecies
        me['_MW'] = np.array([i.MW for i in compounds])
    
    def __setattr__(self, ID, compound):
        if self._immutable:
            raise AttributeError("object is read-only")
        elif isinstance(compound, cmps.Compound):
            compound.ID = ID
            setattr_(self, ID, compound)
        elif ID == 'IDs':
            setattr_(self, '__dict__', {i: getattr(self, i) for i in compound})
        else:
            raise TypeError('can only set Compound objects as attributes')
    
    def __delattr__(self, ID):
        if self._immutable:
            raise AttributeError("object is read-only")
        else:
            super().__delattr__(ID)
    
    def __len__(self):
        return self._Nspecies or len(self.__dict__)
    
    def __iter__(self):
        yield from self._compounds or self.__dict__.values()
    
    @property
    def IDs(self):
        return self._IDs or tuple(self.__dict__.keys())

    def getprops(self, IDs, prop_ID, T, P, phase):
        """Return list of the desired property, prop_ID, for each compound in IDz."""
        props = []; getattr_ = getattr
        for ID in IDs:
            s = getattr_(self, ID)
            s.T = T
            s.P = P
            s.phase = phase
            props.append(getattr_(s, prop_ID))
        return props
    def setprops(self, IDs, prop_ID, new_values):
        """Set new values to compound property, prop_ID, for given IDs."""
        getattr_ = getattr
        for ID in IDs: 
            setattr_(getattr_(self, ID), prop_ID, new_values)

    ### Material Properties ###

    # General methods for getting ideal mixture properties
    def _props(self, ID, mol, T, P, phase):
        """Return component property list."""
        out = np.zeros(self._Nspecies)
        getattr_ = getattr
        ic = self._numcompounds.__iter__().__next__
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
        species = []; index = []; ic = self._numcompounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.UNIFAC_Dortmund_groups and c.Tb not in (inf, None):
                    species.append(c); index.append(i)
            else: ic()
        return species, index

    def _heavy_species(self, mol):
        """Return species and indexes of heavy species not in equilibrium."""
        species = []; index = []; ic = self._numcompounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.Tb in (inf, None) or not c.UNIFAC_Dortmund_groups:
                    species.append(c); index.append(i)
            else: ic()
        return species, index

    def _light_species(self, mol):
        """Return species and indexes of light species not in equilibrium."""
        species = []; index = []; ic = self._numcompounds.__iter__().__next__
        for m in mol:
            if m:
                i, c = ic()
                if c.Tb is -inf:
                    species.append(c); index.append(i)
            else: ic()
        return species, index
    
    # Representation
    def __repr__(self):
        return f'<Species: {", ".join(self.IDs)}>'

    
