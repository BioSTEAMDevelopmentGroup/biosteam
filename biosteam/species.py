# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 13:42:33 2018

@author: yoelr
"""
from biosteam.chemical import Chemical
from biosteam.specie import Specie
from biosteam import np, units_of_measure
from biosteam.exceptions import IDconflict

# TODO: Fix how Dortmund groups are selected in thermo.Chemical. Glycerol defaults to wrong value.

# %% Managing ChEDL thermodynamic properties
inf = np.inf
TPphase = ('T', 'P', 'phase')

class Species:
    """Create Species object that contains biosteam 'Chemical' objects as attributes (by ID), or 'Specie' objects if an ID is invalid.

    **Parameters**

        ***IDs:** [str] Strings should be one of the following [-]:
                  * Name, in IUPAC form or common form or a synonym registered in PubChem
                  * InChI name, prefixed by 'InChI=1S/' or 'InChI=1/'
                  * InChI key, prefixed by 'InChIKey='
                  * PubChem CID, prefixed by 'PubChem='
                  * SMILES (prefix with 'SMILES=' to ensure smiles parsing)
                  * CAS number
         
        **Optional**
        
        **cls:** [type Specie] The type of objects that will be created and stored as attributes. Defaults to biosteam.Chemical
        
    """
    units = units_of_measure

    def __init__(self, *ID, cls=None):
        # For keeping track of IDs/CAS numbers and their order
        super().__setattr__('_ID', [])
        super().__setattr__('_CAS', [])
        
        # Sphinx docs look ugly if a class is in the default argument
        if cls is None:
            cls = Chemical
        
        # Make sure ID is a tuple
        if not ID:
            return
        elif isinstance(ID, str):
            ID = (ID,)
        elif not isinstance(ID[0], str):
            ID = ID[0] # For backwards compatibility
        
        # Set Specie object attributes
        for n in ID:
            n = n.replace('_', ' ')
            try:
                specie = cls(n)
            except:
                raise Exception(f"Chemical, '{n}', not defined in databank.")
            setattr(self, n, specie)

    def __setattr__(self, ID, specie):
        if ID == 'ID':
            super().__setattr__(ID, specie)
        elif isinstance(specie, Specie):
            # Get specie ID and CAS number
            specie.ID = ID = ID.replace(' ', '_')
            CAS = specie.CAS
            
            # Prevent repeated definitions
            IDs = self._ID
            CASs = self._CAS
            if ID in IDs:
                raise IDconflict(f"A Specie object with ID, '{ID}', is already registered.")
            elif CAS in CASs:
                raise IDconflict(f"A Specie object with CAS number, '{CAS}', is already registered.")
            
            # Set attributes and keep track of IDs/CAS numbers
            super().__setattr__(ID, specie)
            super().__setattr__(specie.CAS, specie)
            IDs.append(ID)
            CASs.append(CAS)
        else:
            raise TypeError('Can only set Specie objects as attributes')

    def __delattr__(self, ID):
        # Remove ID and CAS from list
        CAS = getattr(self, ID).CAS
        self._ID.remove(ID)
        self._CAS.remove(CAS)
        
        # Deleat attributes
        super().__delattr__(ID)
        super().__delattr__(CAS)
        
    
    def __iter__(self):
        return (getattr(self, i) for i in self.ID)
    
    @property
    def CAS(self):
        return tuple(self._CAS)
    
    @property
    def ID(self):
        return tuple(self._ID)
    
    @ID.setter
    def ID(self, ID):
        self._ID[:] = ID
        CASs = self._CAS
        CASs.clear()
        for specie in self:
            CASs.append(specie.CAS)
    
    @property
    def _equilibrium_species(self):
        """Return IDs and indeces of species available for equilibrium calculations."""
        specie_IDs = []
        specie_index = []
        
        i = 0
        for sp_ID in self._ID:
            s = getattr(self, sp_ID)
            if s.UNIFAC_groups and s.Tb not in (inf, None):
                specie_IDs.append(sp_ID)
                specie_index.append(i)
            i += 1
        return specie_IDs, specie_index
    
    @property
    def _light_species(self):
        """Return IDs and indeces of species not available for equilibrium calculations."""
        specie_IDs = []
        specie_index = []
        
        i = 0
        for sp_ID in self._ID:
            s = getattr(self, sp_ID)
            if s.Tb is -inf:
                specie_IDs.append(sp_ID)
                specie_index.append(i)
            i += 1
        return specie_IDs, specie_index
    
    @property
    def _heavy_species(self):
        """Return IDs and indeces of species not available for equilibrium calculations."""
        specie_IDs = []
        specie_index = []
        
        i = 0
        for sp_ID in self._ID:
            s = getattr(self, sp_ID)
            if s.Tb in (inf, None) or not s.UNIFAC_groups:
                specie_IDs.append(sp_ID)
                specie_index.append(i)
            i += 1
        return specie_IDs, specie_index
    
    @classmethod
    def combine(cls, *others):
        """Combine Specie\Chemical\Species objects into a new Species object."""
        new = cls()
        for other in others:
            if isinstance(other, Specie):
                setattr(new, other.ID, other)
            elif isinstance(other, Species):
                for i in other.ID:
                    setattr(new, i, getattr(other, i))  
        return new

    def get_props(self, specie_IDs, prop_ID, T=None, P=None, phase=None) -> list:
        """Return list of the desired property, prop_ID, for each specie in species_ID."""
        if isinstance(specie_IDs, str):
            specie = getattr(self, specie_IDs)
            for attr, val in zip(TPphase, (T, P, phase)):
                if val: setattr(specie, attr, val)
            return getattr(specie, prop_ID)
        else:
            props = []
            for ID in specie_IDs:
                specie = getattr(self, ID)
                for attr, val in zip(TPphase, (T, P, phase)):
                    if val: setattr(specie, attr, val)
                props.append(getattr(specie, prop_ID))
        return props

    def set_props(self, specie_IDs, prop_ID, new_values):
        """Set new values to specie property, prop_ID, for each specie in specie IDs."""
        for ID in specie_IDs:
            specie = getattr(self, ID)
            setattr(specie, prop_ID, new_values)

    def _get_Psat(self, specie_IDs, T) -> list:
        """Return list of saturation pressures for the given species."""
        Psat = []
        for ID in specie_IDs:
            Psat.append(getattr(self, ID).VaporPressure(T))
        return Psat

    def _get_Tsat(self, specie_IDs, P) -> list:
        """Return list of temperatures for the given species."""
        Tsat = []
        for ID in specie_IDs:
            Tsat.append(getattr(self, ID).Tsat(P))
        return Tsat

    def _get_Tb(self, specie_IDs) -> list:
        """Return list of normal boiling points for the given species."""
        Tb = []
        for ID in specie_IDs:
            Tb.append(getattr(self, ID).Tb)
        return Tb

    # Representation
    def _info(self):
        """Information on self"""
        return (f'Species:\n {self.ID}').replace("'", '')

    def show(self):
        """print information on self"""
        print(self._info())

    def __repr__(self):
        return (f'<Species: {self.ID}>').replace("'", '')
