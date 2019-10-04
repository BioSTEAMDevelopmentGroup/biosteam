# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

from __future__ import division

__all__ = ['Chemical', 'reference_states']

from math import exp
from .identifiers import CAS_from_any, pubchem_db
from .vapor_pressure import VaporPressure
from .phase_change import Tb, Tm, Hfus, Hsub, EnthalpyVaporization
from .activity import identify_phase
from .critical import Tc, Pc, Vc
from .acentric import omega, StielPolar
from .triple import Tt, Pt
from .thermal_conductivity import ThermalConductivityLiquid, ThermalConductivityGas
from .volume import VolumeGas, VolumeLiquid, VolumeSolid
from .permittivity import Permittivity
from .heat_capacity import HeatCapacitySolid, HeatCapacityGas, HeatCapacityLiquid
from .interface import SurfaceTension
from .viscosity import ViscosityLiquid, ViscosityGas
from .reaction import Hf
from .combustion import Hcombustion
from .safety import Tflash, Tautoignition, LFL, UFL, TWA, STEL, Ceiling, Skin, Carcinogen
from .solubility import SolubilityParameter
from .dipole import dipole_moment as dipole
from .utils import isentropic_exponent, Z, B_from_Z, isobaric_expansion, Joule_Thomson, nu_mu_converter, thermal_diffusivity
from fluids.core import Reynolds, Capillary, Weber, Bond, Grashof, Peclet_heat
from .lennard_jones import Stockmayer, MolecularDiameter
from .environment import GWP, ODP, logP
from .law import legal_status, economic_status
from .refractivity import refractive_index
from .electrochem import conductivity
from .elements import atom_fractions, mass_fractions, similarity_variable, atoms_to_Hill, simple_formula_parser, molecular_weight, charge_from_formula
from .eos import GCEOS_DUMMY, PR
from .unifac import DDBST_UNIFAC_assignments, DDBST_MODIFIED_UNIFAC_assignments, DDBST_PSRK_assignments, load_group_assignments_DDBST, UNIFAC_RQ, Van_der_Waals_volume, Van_der_Waals_area
from fluids.core import Prandtl, Parachor, Jakob
from scipy.optimize import newton
from scipy.constants import physical_constants
from math import log

_R = 8.3144598  # Universal gas constant (J/mol-K)

# RDKIT
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    # pragma: no cover
    pass


#import warnings
#warnings.filterwarnings("ignore")

caching = True

# Format: (T, P, phase, H, S, molar=True)
IAPWS = (273.16, 611.655, 'l', 0.00922, 0, True) # Water; had to convert Href from mass to molar
ASHRAE = (233.15, 'Psat', 'l', 0, 0, True) # As described in REFPROP
IIR = (273.15, 'Psat', 'l', 200E3, 1000, False) # 200 kj/kg reference, as described in REFPROP
REFPROP = ('Tb', 101325, 'l', 0, 0, True)
CHEMSEP = (298., 101325, 'g', 0, 0, True) # It has an option to add Hf to the reference
PRO_II = (298.15, 101325, 'gas', 0, 0, True)
HYSYS = (298.15, 101325, 'calc', 'Hf', 0, True)
UNISIM = HYSYS #
SUPERPRO = (298.15, 101325, 'calc', 0, 0, True) # No support for entropy found, 0 assumed
# note soecifying a phase works for chemicals but not mixtures.

reference_states = [IAPWS, ASHRAE, IIR, REFPROP, CHEMSEP, PRO_II, HYSYS,
                    UNISIM, SUPERPRO]
_chemical_cache = {}

def phase_property(func):
    attr = func.__name__
    phaseprop = lambda self: getattr(self, attr + self.phase)
    phaseprop.__name__ = attr
    phaseprop.__doc__ = func.__doc__
    return property(phaseprop)


class Chemical:
    '''Creates a Chemical object which contains basic information such as 
    molecular weight and the structure of the species, as well as thermodynamic
    and transport properties as a function of temperature and pressure.
    
    Parameters
    ----------
    ID : str
        One of the following [-]:
            * Name, in IUPAC form or common form or a synonym registered in PubChem
            * InChI name, prefixed by 'InChI=1S/' or 'InChI=1/'
            * InChI key, prefixed by 'InChIKey='
            * PubChem CID, prefixed by 'PubChem='
            * SMILES (prefix with 'SMILES=' to ensure smiles parsing)
            * CAS number
    T : float, optional
        Temperature of the chemical (default 298.15 K), [K]
    P : float, optional
        Pressure of the chemical (default 101325 Pa) [Pa]
    
    Attributes
    ----------
    T : float
        Temperature of the chemical, [K]
    P : float
        Pressure of the chemical, [Pa]
    phase : str
        Phase of the chemical; one of 's', 'l', 'g', or 'l/g'.
    ID : str
        User specified string by which the chemical's CAS was looked up.
    CAS : str
        The CAS number of the chemical.
    PubChem : int
        PubChem Compound identifier (CID) of the chemical; all chemicals are
        sourced from their database. Chemicals can be looked at online at
        `<https://pubchem.ncbi.nlm.nih.gov>`_.
    MW : float
        Molecular weight of the compound, [g/mol]
    formula : str
        Molecular formula of the compound.
    atoms : dict
        dictionary of counts of individual atoms, indexed by symbol with
        proper capitalization, [-]
    similarity_variable : float
        Similarity variable, see :obj:`thermo.elements.similarity_variable`
        for the definition, [mol/g]
    smiles : str
        Simplified molecular-input line-entry system representation of the
        compound.
    InChI : str
        IUPAC International Chemical Identifier of the compound.
    InChI_Key : str
        25-character hash of the compound's InChI.
    IUPAC_name : str
        Preferred IUPAC name for a compound.
    synonyms : list of strings
        All synonyms for the compound found in PubChem, sorted by popularity.
    Tm : float
        Melting temperature [K]
    Tb : float
        Boiling temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure [Pa]
    Vc : float
        Critical volume [m^3/mol]
    Zc : float
        Critical compressibility [-]
    rhoc : float
        Critical molar density [mol/m^3]
    omega : float
        Acentric factor [-]
    StielPolar : float
        Stiel Polar factor, see :obj:`thermo.acentric.StielPolar` for
        the definition [-]
    Tt : float
        Triple temperature, [K]
    Pt : float
        Triple pressure, [Pa]
    Hfus : float
        Molar enthalpy of fusion [J/mol]
    Hsub : float
        Molar enthalpy of sublimation [J/mol]
    Hf : float
        Enthalpy of formation [J/mol]
    Hc : float
        Molar enthalpy of combustion [J/mol]
    Tflash : float
        Flash point of the chemical, [K]
    Tautoignition : float
        Autoignition point of the chemical, [K]
    LFL : float
        Lower flammability limit of the gas in an atmosphere at STP, mole 
        fraction [-]
    UFL : float
        Upper flammability limit of the gas in an atmosphere at STP, mole 
        fraction [-]
    TWA : tuple[quantity, unit]
        Time-Weighted Average limit on worker exposure to dangerous chemicals.
    STEL : tuple[quantity, unit]
        Short-term Exposure limit on worker exposure to dangerous chemicals.
    Ceiling : tuple[quantity, unit]
        Ceiling limits on worker exposure to dangerous chemicals.
    Skin : bool
        Whether or not a chemical can be absorbed through the skin.
    Carcinogen : str or dict
        Carcinogen status information.
    dipole : float
        Dipole moment in debye, [3.33564095198e-30 ampere*second^2]
    Stockmayer : float
        Lennard-Jones depth of potential-energy minimum over k, [K]
    MolecularDiameter : float
        Lennard-Jones molecular diameter, [angstrom]
    GWP : float
        Global warming potential (default 100-year outlook) (impact/mass 
        chemical)/(impact/mass CO2), [-]
    ODP : float
        Ozone Depletion potential (impact/mass chemical)/(impact/mass CFC-11),
        [-]
    logP : float
        Octanol-water partition coefficient, [-]
    legal_status : str or dict
        Legal status information [-]
    economic_status : list
        Economic status information [-]
    RI : float
        Refractive Index on the Na D line, [-]
    RIT : float
        Temperature at which refractive index reading was made
    conductivity : float
        Electrical conductivity of the fluid, [S/m]
    conductivityT : float
        Temperature at which conductivity measurement was made
    VaporPressure : TDependentModel
    EnthalpyVaporization : TPDependentModel
    VolumeSolid : TPDependentModel
    VolumeLiquid : TPDependentModel
    VolumeGas : TPDependentModel
    HeatCapacitySolid : TDependentModel
    HeatCapacityLiquid : TDependentModel
    HeatCapacityGas : TDependentModel
    ViscosityLiquid : TDependentModel
    ViscosityGas : TDependentModel
    ThermalConductivityLiquid : TDependentModel
    ThermalConductivityGas : TDependentModel
    SurfaceTension : TDependentModel
    Permittivity : TDependentModel
    SolubilityParameter : TDependentModel
    Psat_298 : float
        Vapor pressure of the chemical at 298.15 K, [Pa]
    phase_STP : str
        Phase of the chemical at 298.15 K and 101325 Pa; one of 's', 'l', 'g',
        or 'l/g'.
    Vl_Tb : float
        Molar volume of liquid phase at the normal boiling point [m^3/mol]
    Vl_Tm : float
        Molar volume of liquid phase at the melting point [m^3/mol]
    Vl_STP : float
        Molar volume of liquid phase at 298.15 K and 101325 Pa [m^3/mol]
    Vg_STP : float
        Molar volume of gas phase at 298.15 K and 101325 Pa [m^3/mol]
    Hvap_Tb : float
        Molar enthalpy of vaporization at the normal boiling point [J/mol]
    alpha
    alphag
    alphal
    API
    aromatic_rings
    atom_fractions
    Bvirial
    charge
    Cp
    Cpg
    Cpg
    Cpl
    Cpl
    Cp
    Cps
    Cvg
    eos
    Hill
    Hvap
    isentropic_exponent
    isobaric_expansion
    isobaric_expansion_g
    isobaric_expansion_l
    JT
    JTg
    JTl
    k
    kg
    kl
    mass_fractions
    mu
    mug
    mul
    nu
    nug
    nul
    Parachor
    permittivity
    Poynting
    Pr
    Prg
    Prl
    Psat
    PSRK_groups
    rdkitmol
    rdkitmol_Hs
    rho
    rhog
    rhol
    rhos
    rings
    SG
    SGg
    SGl
    SGs
    sigma
    UNIFAC_Dortmund_groups
    UNIFAC_groups
    UNIFAC_R
    UNIFAC_Q
    Van_der_Waals_area
    Van_der_Waals_volume
    V
    Vg
    Vl
    Vs
    Z
    Zg
    Zl
    Zs
    '''

    __atom_fractions = None
    __mass_fractions = None
    __UNIFAC_groups = None
    __UNIFAC_Dortmund_groups = None
    __PSRK_groups = None
    __rdkitmol = None
    __rdkitmol_Hs = None
    __Hill = None
    __legal_status = None
    __economic_status = None
    def __repr__(self):
        return '<Chemical [%s], T=%.2f K, P=%.0f Pa>' %(self.name, self.T, self.P)

    def __init__(self, ID, T=298.15, P=101325):
        if isinstance(ID, dict):
            self.CAS = ID['CASRN']
            self.ID = self.name = ID['name']
            self.formula = ID['formula']
            self.MW = ID['MW'] if 'MW' in ID else molecular_weight(simple_formula_parser(self.formula))
            self.PubChem = ID['PubChem'] if 'PubChem' in ID else None
            self.smiles = ID['smiles'] if 'smiles' in ID else None
            self.InChI = ID['InChI'] if 'InChI' in ID else None
            self.InChI_Key = ID['InChI_Key'] if 'InChI_Key' in ID else None
            self.synonyms = ID['synonyms'] if 'synonyms' in ID else None
        else:
            self.ID = ID
            # Identification
            self.CAS = CAS_from_any(ID)
            self.ChemicalMetadata = pubchem_db.search_CAS(self.CAS)
        self.T = T
        self.P = P
        if self.CAS in _chemical_cache and caching:
            self.__dict__.update(_chemical_cache[self.CAS].__dict__)
            self.refresh()
        else:
            if not isinstance(ID, dict):
                self.PubChem = self.ChemicalMetadata.pubchemid
                self.MW = self.ChemicalMetadata.MW
                self.formula = self.ChemicalMetadata.formula
                self.smiles = self.ChemicalMetadata.smiles
                self.InChI = self.ChemicalMetadata.InChI
                self.InChI_Key = self.ChemicalMetadata.InChI_key
                self.IUPAC_name = self.ChemicalMetadata.iupac_name.lower()
                self.name = self.ChemicalMetadata.common_name.lower()
                self.synonyms = self.ChemicalMetadata.all_names
            self.atoms = simple_formula_parser(self.formula)
            self.similarity_variable = similarity_variable(self.atoms, self.MW)
            self.eos_in_a_box = []
            self.set_constants()
            self.set_eos(T=T, P=P)
            self.set_models()
            self.set_ref()
            self.refresh()
            if len(_chemical_cache) < 1000:
                _chemical_cache[self.CAS] = self

    def refresh(self):
        self.phase = identify_phase(T=self.T, P=self.P, Tm=self.Tm, Tb=self.Tb, Tc=self.Tc, Psat=self.Psat)
        self.eos = self.eos.to_TP(T=self.T, P=self.P)
        self.eos_in_a_box[0] = self.eos
        self.set_thermo()

    def draw_2d(self, width=300, height=300, Hs=False): # pragma: no cover
        r'''Interface for drawing a 2D image of the molecule.
        Requires an HTML5 browser, and the libraries RDKit and
        IPython. An exception is raised if either of these libraries is
        absent.

        Parameters
        ----------
        width : int
            Number of pixels wide for the view
        height : int
            Number of pixels tall for the view
        Hs : bool
            Whether or not to show hydrogen

        Examples
        --------
        >>> Chemical('decane').draw_2d() # doctest: +ELLIPSIS
        <PIL.Image.Image image mode=RGBA size=300x300 at 0x...>
        '''
        try:
            from rdkit.Chem import Draw
            if Hs:
                mol = self.rdkitmol_Hs
            else:
                mol = self.rdkitmol
            return Draw.MolToImage(mol, size=(width, height))
        except:
            return 'Rdkit is required for this feature.'

    def draw_3d(self, width=300, height=500, style='stick', Hs=True): # pragma: no cover
        r'''Interface for drawing an interactive 3D view of the molecule.
        Requires an HTML5 browser, and the libraries RDKit, pymol3D, and
        IPython. An exception is raised if all three of these libraries are
        not installed.

        Parameters
        ----------
        width : int
            Number of pixels wide for the view
        height : int
            Number of pixels tall for the view
        style : str
            One of 'stick', 'line', 'cross', or 'sphere'
        Hs : bool
            Whether or not to show hydrogen

        Examples
        --------
        >>> Chemical('cubane').draw_3d()
        <IPython.core.display.HTML object>
        '''
        try:
            import py3Dmol
            from IPython.display import display
            if Hs:
                mol = self.rdkitmol_Hs
            else:
                mol = self.rdkitmol
            AllChem.EmbedMultipleConfs(mol)
            mb = Chem.MolToMolBlock(mol)
            p = py3Dmol.view(width=width,height=height)
            p.addModel(mb,'sdf')
            p.setStyle({style:{}})
            p.zoomTo()
            display(p.show())
        except:
            return 'py3Dmol, RDKit, and IPython are required for this feature.'

    def set_constants(self):
        CAS = self.CAS
        MW = self.MW
        self.Tm = Tm(CAS)
        self.Tb = Tb(CAS)

        # Critical Point
        self.Tc = Tc(CAS)
        self.Pc = Pc(CAS)
        self.Vc = Vc(CAS)
        self.omega = omega(CAS)
        self.StielPolar = StielPolar(CAS, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        self.Zc = Z(self.Tc, self.Pc, self.Vc) if all((self.Tc, self.Pc, self.Vc)) else None
        self.rhoc = 1./self.Vc if self.Vc else None

        # Triple point
        self.Pt = Pt(CAS)
        self.Tt = Tt(CAS)

        # Enthalpy
        self.Hfus = Hfus(CAS, MW=MW)
        self.Hsub = Hsub(CAS, MW=MW)

        # Chemistry
        self.Hf = Hf(CAS)
        self.Hc = Hcombustion(atoms=self.atoms, Hf=self.Hf)

        # Fire Safety Limits
        self.Tflash = Tflash(CAS)
        self.Tautoignition = Tautoignition(CAS)
        self.LFL = LFL(CAS, atoms=self.atoms, Hc=self.Hc)
        self.UFL = UFL(CAS, atoms=self.atoms, Hc=self.Hc)

        # Chemical Exposure Limits
        self.TWA = TWA(CAS)
        self.STEL = STEL(CAS)
        self.Ceiling = Ceiling(CAS)
        self.Skin = Skin(CAS)
        self.Carcinogen = Carcinogen(CAS)

        # Misc
        self.dipole = dipole(CAS) # Units of Debye
        self.Stockmayer = Stockmayer(CAS, Tm=self.Tm, Tb=self.Tb, Tc=self.Tc,
                                     Zc=self.Zc, omega=self.omega)

        # Environmental
        self.GWP = GWP(CAS)
        self.ODP = ODP(CAS)
        self.logP = logP(CAS)

        # Analytical
        self.RI, self.RIT = refractive_index(CAS)
        self.conductivity, self.conductivityT = conductivity(CAS)

    def set_eos(self, T, P, eos=PR):
        try:
            self.eos = eos(T=T, P=P, Tc=self.Tc, Pc=self.Pc, omega=self.omega)
        except:
            # Handle overflow errors and so on
            self.eos = GCEOS_DUMMY(T=T, P=P)

    @property
    def eos(self):
        r'''Equation of state object held by the chemical; used to calculate
        excess thermodynamic quantities, and also provides a vapor pressure
        curve, enthalpy of vaporization curve, fugacity, thermodynamic partial
        derivatives, and more; see :obj:`thermo.eos` for a full listing.

        Examples
        --------
        >>> Chemical('methane').eos.V_g
        0.02441019502181826
        '''
        return self.eos_in_a_box[0]

    @eos.setter
    def eos(self, eos):
        if self.eos_in_a_box: self.eos_in_a_box.pop()
        # Pass this mutable list to objects so if it is changed, it gets
        # changed in the property method too
        self.eos_in_a_box.append(eos)

    def set_models(self):
        # Tempearture and Pressure Denepdence
        # Get and choose initial methods
        eos = self.eos_in_a_box
        CAS = self.CAS
        MW = self.MW
        Tm = self.Tm
        Tb = self.Tb
        Tc = self.Tc
        Pc = self.Pc
        Vc = self.Vc
        Zc = self.Zc
        omega = self.omega
        dipole = self.dipole
        similarity_variable = self.similarity_variable
        self.VaporPressure = VaporPressure(CAS, Tb=Tb, Tc=Tc, Pc=Pc,
                                           omega=omega,
                                           eos=eos)
        Psat = self.VaporPressure.evaluate
        self.Psat_298 = Psat(298.15)
        self.phase_STP = identify_phase(T=298.15, P=101325., Tm=Tm,
                                        Tb=Tb, Tc=Tc, Psat=self.Psat_298)

        self.VolumeLiquid = VolumeLiquid(CAS, MW=MW, Tb=Tb,
                                         Tc=Tc, Pc=Pc, Vc=Vc, Zc=Zc,
                                         omega=omega,
                                         dipole=self.dipole,
                                         Psat=Psat,
                                         eos=eos)

        self.Vl_Tb = self.VolumeLiquid.evaluate(Tb) if Tb else None
        self.Vl_Tm = self.VolumeLiquid.evaluate(Tm) if Tm else None
        self.Vl_STP = self.VolumeLiquid.evaluate(298.15)

        self.VolumeGas = VolumeGas(CAS, MW=MW, Tc=Tc, Pc=Pc,
                                   omega=omega, dipole=self.dipole,
                                   eos=self.eos_in_a_box)
        self.Vg_STP = self.VolumeGas.evaluate(298.15, 101325)
        self.VolumeSolid = VolumeSolid(CAS, MW=MW, Tt=self.Tt)
        self.HeatCapacityGas = HeatCapacityGas(CAS, MW=MW,
                                               similarity_variable=similarity_variable)

        self.HeatCapacitySolid = HeatCapacitySolid(CAS, MW=MW,
                                                   similarity_variable=similarity_variable)

        self.HeatCapacityLiquid = HeatCapacityLiquid(CAS, MW=self.MW,
                                                     similarity_variable=self.similarity_variable,
                                                     Cpgm=self.HeatCapacityGas.evaluate,
                                                     Tc=Tc, omega=omega)

        self.EnthalpyVaporization = EnthalpyVaporization(CAS, Tb=Tb, Tc=Tc,
                                                         Pc=Pc, omega=omega,
                                                         similarity_variable=similarity_variable)
        self.Hvap_Tb = self.EnthalpyVaporization.evaluate(Tb) if Tb else None

        self.ViscosityLiquid = ViscosityLiquid(CAS, MW=MW, Tm=Tm, Tc=Tc,
                                               Pc=Pc, Vc=Vc, omega=omega,
                                               Psat=Psat,
                                               Vl=self.VolumeLiquid.evaluate)

        Vg_atm_T_dependent = lambda T : self.VolumeGas.evalute(T, 101325)
        self.ViscosityGas = ViscosityGas(CAS, MW=MW, Tc=Tc, Pc=Pc,
                                         Zc=Zc, dipole=dipole, Vg=Vg_atm_T_dependent)

        self.ThermalConductivityLiquid = ThermalConductivityLiquid(CAS, MW=MW, Tm=Tm, Tb=Tb, Tc=Tc,
                                                                   Pc=Pc, omega=omega, Hfus=self.Hfus)

        Cvgm_calc = lambda T : self.HeatCapacityGas.evaluate(T) - _R
        self.ThermalConductivityGas = ThermalConductivityGas(CAS, MW=MW, Tb=Tb, Tc=Tc,
                                                             Pc=Pc, Vc=Vc, Zc=Zc, omega=omega,
                                                             dipole=self.dipole, Vg=self.VolumeGas,
                                                             Cvgm=Cvgm_calc, mug=self.ViscosityGas.evaluate)
        self.SurfaceTension = SurfaceTension(CAS, MW=self.MW, Tb=self.Tb,
                                             Tc=self.Tc, Pc=self.Pc, Vc=self.Vc, Zc=self.Zc,
                                             omega=self.omega, StielPolar=self.StielPolar,
                                             Hvap_Tb=self.Hvap_Tb, Vl=self.VolumeLiquid.evaluate,
                                             Cpl=self.HeatCapacityLiquid.evaluate)
        self.Permittivity = Permittivity(CAS)
        self.SolubilityParameter = SolubilityParameter(CAS, Hvapm=self.Hvap_Tbm, Vl=self.Vl_STP,)
        self.MolecularDiameter = MolecularDiameter(CAS, Tc=Tc, Pc=Pc, Vc=Vc, Zc=Zc, omega=omega,
                                                   V=self.Vl_Tm, Vb=self.Vl_Tb)


    def set_ref(self, T_ref=298.15, P_ref=101325, phase_ref='calc', H_ref=0, S_ref=0):
        # Muse run after set_models, set_phase due to HeatCapacity*, phase_STP
        self.T_ref = getattr(self, T_ref) if isinstance(T_ref, str) else T_ref
        self.P_ref = getattr(self, P_ref) if isinstance(P_ref, str) else P_ref
        self.H_ref = getattr(self, H_ref) if isinstance(H_ref, str) else H_ref
        self.S_ref = getattr(self, S_ref) if isinstance(S_ref, str) else S_ref
        self.phase_ref = self.phase_STP if phase_ref == 'calc' else phase_ref
        Cps_int = self.HeatCapacitySolid.integrate
        Cpl_int = self.HeatCapacityLiquid.integrate
        Cpg_int = self.HeatCapacityGas.integrate
        Cps_int_T = self.HeatCapacitySolid.integrate_over_T
        Cpl_int_T = self.HeatCapacityLiquid.integrate_over_T
        Cpg_int_T = self.HeatCapacityGas.integrate_over_T

        # Integrals stored to avoid recalculation, all from T_low to T_high
        try:
            # Enthalpy integrals
            if self.phase_ref != 'l' and self.Tm and self.Tb:
                self.H_int_l_Tm_to_Tb = Cpl_int(self.Tm, self.Tb)
            if self.phase_ref == 's' and self.Tm:
                self.H_int_T_ref_s_to_Tm = Cps_int(self.T_ref, self.Tm)
            if self.phase_ref == 'g' and self.Tb:
                self.H_int_Tb_to_T_ref_g = Cpg_int(self.Tb, self.T_ref)
            if self.phase_ref == 'l' and self.Tm and self.Tb:
                self.H_int_l_T_ref_l_to_Tb = Cpl_int(self.T_ref, self.Tb)
                self.H_int_l_Tm_to_T_ref_l = Cpl_int(self.Tm, self.T_ref)

            # Entropy integrals
            if self.phase_ref != 'l' and self.Tm and self.Tb:
                self.S_int_l_Tm_to_Tb = Cpl_int_T(self.Tm, self.Tb)
            if self.phase_ref == 's' and self.Tm:
                self.S_int_T_ref_s_to_Tm = Cps_int_T(self.T_ref, self.Tm)
            if self.phase_ref == 'g' and self.Tb:
                self.S_int_Tb_to_T_ref_g = Cpg_int_T(self.Tb, self.T_ref)
            if self.phase_ref == 'l' and self.Tm and self.Tb:
                self.S_int_l_T_ref_l_to_Tb = Cpl_int_T(self.T_ref, self.Tb)
                self.S_int_l_Tm_to_T_ref_l = Cpl_int_T(self.Tm, self.T_ref)
        except:
            pass

        # Excess properties stored
        try:
            if self.phase_ref == 'g':
                self.eos_phase_ref = self.eos.to_TP(self.T_ref, self.P_ref)
                self.H_dep_ref_g = self.eos_phase_ref.H_dep_g
                self.S_dep_ref_g = self.eos_phase_ref.S_dep_g


            elif self.phase_ref == 'l':
                self.eos_phase_ref = self.eos.to_TP(self.T_ref, self.P_ref)
                self.H_dep_ref_l = self.eos_phase_ref.H_dep_l
                self.S_dep_ref_l = self.eos_phase_ref.S_dep_l

                self.H_dep_T_ref_Pb = self.eos.to_TP(self.T_ref, 101325).H_dep_l
                self.S_dep_T_ref_Pb = self.eos.to_TP(self.T_ref, 101325).S_dep_l

            if self.Tb:
                self.eos_Tb = self.eos.to_TP(self.Tb, 101325)
                self.H_dep_Tb_Pb_g = self.eos_Tb.H_dep_g
                self.H_dep_Tb_Pb_l = self.eos_Tb.H_dep_l

                self.H_dep_Tb_P_ref_g = self.eos.to_TP(self.Tb, self.P_ref).H_dep_g
                self.S_dep_Tb_P_ref_g = self.eos.to_TP(self.Tb, self.P_ref).S_dep_g


                self.S_dep_Tb_Pb_g = self.eos_Tb.S_dep_g
                self.S_dep_Tb_Pb_l = self.eos_Tb.S_dep_l
        except:
            pass

    @property
    def H_excess(self):
        H_dep = 0
        T = self.T
        P = self.P
        if self.phase_ref == 'g' and self.phase == 'g':
            H_dep += self.eos.to_TP(T, P).H_dep_g - self.H_dep_ref_g

        elif self.phase_ref == 'l' and self.phase == 'l':
            try:
                H_dep += self.eos.to_TP(T, P).H_dep_l - self._eos_T_101325.H_dep_l
            except:
                H_dep += 0

        elif self.phase_ref == 'g' and self.phase == 'l':
            H_dep += self.H_dep_Tb_Pb_g - self.H_dep_Tb_P_ref_g
            H_dep += (self.eos.to_TP(T, P).H_dep_l - self._eos_T_101325.H_dep_l)

        elif self.phase_ref == 'l' and self.phase == 'g':
            H_dep += self.H_dep_T_ref_Pb - self.H_dep_ref_l
            H_dep += (self.eos.to_TP(T, P).H_dep_g - self.H_dep_Tb_Pb_g)
        return H_dep

    @property
    def S_excess(self):
        S_dep = 0
        T = self.T
        P = self.P
        if self.phase_ref == 'g' and self.phase == 'g':
            S_dep += self.eos.to_TP(T, P).S_dep_g - self.S_dep_ref_g

        elif self.phase_ref == 'l' and self.phase == 'l':
            try:
                S_dep += self.eos.to_TP(T, P).S_dep_l - self._eos_T_101325.S_dep_l
            except:
                S_dep += 0

        elif self.phase_ref == 'g' and self.phase == 'l':
            S_dep += self.S_dep_Tb_Pb_g - self.S_dep_Tb_P_ref_g
            S_dep += (self.eos.to_TP(T, P).S_dep_l - self._eos_T_101325.S_dep_l)

        elif self.phase_ref == 'l' and self.phase == 'g':
            S_dep += self.S_dep_T_ref_Pb - self.S_dep_ref_l
            S_dep += (self.eos.to_TP(T, P).S_dep_g - self.S_dep_Tb_Pb_g)
        return S_dep

    ### Energies ###

    @property
    def H(self):
        """Enthapy (kJ/kmol) disregarding pressure and assuming the specified phase."""
        phase_ref = self.phase_ref
        phase = self.phase

        # Perform enthalpy calculations at 101325 Pa
        if phase_ref == phase:
            if phase == 'l':
                H = self.HeatCapacityLiquid.integrate(self.T_ref, self.T)
            elif phase == 'g':
                H = self.HeatCapacityGas.integrate(self.T_ref, self.T)
            elif phase == 's':
                H = self.HeatCapacitySolid.integrate(self.T_ref, self.T)
        elif phase_ref == 'l' and phase == 'g':
            H = self.H_int_l_T_ref_l_to_Tb + self.Hvap_Tb + \
                self.HeatCapacityGas.integrate(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 'l':
            H = -self.H_int_Tb_to_T_ref_g - self.Hvap_Tb + \
                self.HeatCapacityLiquid.integrate(self.Tb, self.T)
        elif phase_ref == 's' and phase == 'l':
            H = self.H_int_T_ref_s_to_Tm + self.Hfus + \
                self.HeatCapacityLiquid.integrate(self.Tm, self.T)
        elif phase_ref == 'l' and phase == 's':
            H = -self.H_int_l_Tm_to_T_ref_l - self.Hfus + \
                self.HeatCapacitySolid.integrate(self.Tm, self.T)
        elif phase_ref == 's' and phase == 'g':
            H = self.H_int_T_ref_s_to_Tm + self.Hfus + self.H_int_l_Tm_to_Tb + \
                self.Hvap_Tb + self.HeatCapacityGas.integrate(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 's':
            H = -self.H_int_Tb_to_T_ref_g - self.Hvap_Tb - self.H_int_l_Tm_to_Tb - \
                self.Hfus + self.HeatCapacitySolid.integrate(self.Tm, self.T)
        return self.H_ref + H

    @property
    def S(self):
        """Entropy (kJ/kmol) assuming the specified phase."""
        S = self.S_ref
        phase = self.phase
        phase_ref = self.phase_ref

        # Add Pressure Entropy
        if self.phase == 'g':
            S += -_R*log(self.P/self.P_ref)

        # Add Temperature Entropy
        if phase == phase_ref:
            if phase == 'l':
                S += self.HeatCapacityLiquid.integrate_over_T(self.T_ref, self.T)
            elif phase == 'g':
                S += self.HeatCapacityGas.integrate_over_T(self.T_ref, self.T)
            elif phase == 's':
                S += self.HeatCapacitySolid.integrate_over_T(self.T_ref, self.T)
        elif phase_ref == 'l' and phase == 'g':
            S += self.S_int_l_T_ref_l_to_Tb + self.Hvap_Tb / \
                self.Tb + self.HeatCapacityGas.integrate_over_T(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 'l':
            S += - self.S_int_Tb_to_T_ref_g - self.Hvap / \
                self.Tb + self.HeatCapacityLiquid.integrate_over_T(self.Tb, self.T)
        elif phase_ref == 's' and phase == 'l':
            S += self.S_int_T_ref_s_to_Tm + self.Hfus / \
                self.Tm + self.HeatCapacityLiquid.integrate_over_T(self.Tm, self.T)
        elif self.phase_ref == 'l' and phase == 's':
            S += - self.S_int_l_Tm_to_T_ref_l - self.Hfus / \
                self.Tm + self.HeatCapacitySolid.integrate_over_T(self.Tm, self.T)
        elif phase_ref == 's' and phase == 'g':
            S += self.S_int_T_ref_s_to_Tm + self.Hfus/self.Tm + \
                self.S_int_l_Tm_to_Tb + self.Hvap_Tb / \
                self.Tb + self.HeatCapacityGas.integrate_over_T(self.Tb, self.T)
        elif phase_ref == 'g' and phase == 's':
            S += - self.S_int_Tb_to_T_ref_g - self.Hvap_Tb/self.Tb - \
                self.S_int_l_Tm_to_Tb - self.Hfus / \
                self.Tm + self.HeatCapacitySolid.integrate_over_T(self.Tm, self.T)
        else:
            raise Exception(f'Error in Compound object "{self.ID}" with phase "{phase}" and reference "{phase_ref}.')
        return S

    @property
    def U(self):
        r'''Internal energy of the chemical at its current temperature and
        pressure, in units of [J/mol].

        This property requires that :obj:`thermo.chemical.set_thermo` ran
        successfully to be accurate.
        It also depends on the molar volume of the chemical at its current
        conditions.
        '''
        return self.H - self.P*self.V
    
    @property
    def A(self):
        r'''Helmholtz energy of the chemical at its current temperature and
        pressure, in units of [J/mol].

        This property requires that :obj:`thermo.chemical.set_thermo` ran
        successfully to be accurate.
        It also depends on the molar volume of the chemical at its current
        conditions.
        '''
        return self.U - self.T*self.S

    @property
    def G(self):
        """Gibbs free energy flow rate as a function of T and P, excluding formation energies [kJ/hr].
        
        This property requires that :obj:`thermo.chemical.set_thermo` ran
        successfully to be accurate.
        It also depends on the molar volume of the chemical at its current
        conditions.
        """
        return self.H - self.S*self.T

    def set_TH(self, T, H):
        self.T = T
        def to_solve(P):
            self.P = P
            self.calculate()
            return self.H - H
        newton(to_solve, self.P)

    def set_PH(self, P, H):
        self.P = P
        def to_solve(T):
            self.calculate(T, P)
            return self.H - H
        newton(to_solve, self.T)

    def set_TS(self, T, S):
        self.T = T
        def to_solve(P):
            self.P = P
            self.calculate()
            return self.S - S
        newton(to_solve, self.P)

    def set_PS(self, P, S):
        self.P = P
        def to_solve(T):
            self.calculate(T, P)
            return self.S - S
        newton(to_solve, self.T)

    def set_thermo(self):
        self._eos_T_101325 = self.eos.to_TP(self.T, 101325)

    ### Temperature independent properties - calculate lazily
    @property
    def charge(self):
        r'''Charge of a chemical, computed with RDKit from a chemical's SMILES.
        If RDKit is not available, holds None.

        Examples
        --------
        >>> Chemical('sodium ion').charge
        1
        '''
        try:
            if not self.rdkitmol:
                return charge_from_formula(self.formula)
            else:
                return Chem.GetFormalCharge(self.rdkitmol)
        except:
            return charge_from_formula(self.formula)

    @property
    def rings(self):
        r'''Number of rings in a chemical, computed with RDKit from a
        chemical's SMILES. If RDKit is not available, holds None.

        Examples
        --------
        >>> Chemical('Paclitaxel').rings
        7
        '''
        try:
            return Chem.Descriptors.RingCount(self.rdkitmol)
        except:
            return None

    @property
    def aromatic_rings(self):
        r'''Number of aromatic rings in a chemical, computed with RDKit from a
        chemical's SMILES. If RDKit is not available, holds None.

        Examples
        --------
        >>> Chemical('Paclitaxel').aromatic_rings
        3
        '''
        try:
            return Chem.Descriptors.NumAromaticRings(self.rdkitmol)
        except:
            return None

    @property
    def rdkitmol(self):
        r'''RDKit object of the chemical, without hydrogen. If RDKit is not
        available, holds None.

        For examples of what can be done with RDKit, see
        `their website <http://www.rdkit.org/docs/GettingStartedInPython.html>`_.
        '''
        if self.__rdkitmol:
            return self.__rdkitmol
        else:
            try:
                self.__rdkitmol = Chem.MolFromSmiles(self.smiles)
                return self.__rdkitmol
            except:
                return None

    @property
    def rdkitmol_Hs(self):
        r'''RDKit object of the chemical, with hydrogen. If RDKit is not
        available, holds None.

        For examples of what can be done with RDKit, see
        `their website <http://www.rdkit.org/docs/GettingStartedInPython.html>`_.
        '''
        if self.__rdkitmol_Hs:
            return self.__rdkitmol_Hs
        else:
            try:
                self.__rdkitmol_Hs = Chem.AddHs(self.rdkitmol)
                return self.__rdkitmol_Hs
            except:
                return None

    @property
    def Hill(self):
        r'''Hill formula of a compound. For a description of the Hill system,
        see :obj:`thermo.elements.atoms_to_Hill`.

        Examples
        --------
        >>> Chemical('furfuryl alcohol').Hill
        'C5H6O2'
        '''
        if self.__Hill:
            return self.__Hill
        else:
            self.__Hill = atoms_to_Hill(self.atoms)
            return self.__Hill

    @property
    def atom_fractions(self):
        r'''Dictionary of atom:fractional occurence of the elements in a
        chemical. Useful when performing element balances. For mass-fraction
        occurences, see :obj:`mass_fractions`.

        Examples
        --------
        >>> Chemical('Ammonium aluminium sulfate').atom_fractions
        {'H': 0.25, 'S': 0.125, 'Al': 0.0625, 'O': 0.5, 'N': 0.0625}
        '''
        if self.__atom_fractions:
            return self.__atom_fractions
        else:
            self.__atom_fractions = atom_fractions(self.atoms)
            return self.__atom_fractions

    @property
    def mass_fractions(self):
        r'''Dictionary of atom:mass-weighted fractional occurence of elements.
        Useful when performing mass balances. For atom-fraction occurences, see
        :obj:`atom_fractions`.

        Examples
        --------
        >>> Chemical('water').mass_fractions
        {'H': 0.11189834407236524, 'O': 0.8881016559276347}
        '''
        if self.__mass_fractions:
            return self.__mass_fractions
        else:
            self.__mass_fractions =  mass_fractions(self.atoms, self.MW)
            return self.__mass_fractions

    @property
    def legal_status(self):
        r'''Dictionary of legal status indicators for the chemical.

        Examples
        --------
        >>> pprint(Chemical('benzene').legal_status)
        {'DSL': 'LISTED',
         'EINECS': 'LISTED',
         'NLP': 'UNLISTED',
         'SPIN': 'LISTED',
         'TSCA': 'LISTED'}
        '''
        if self.__legal_status:
            return self.__legal_status
        else:
            self.__legal_status = legal_status(self.CAS, Method='COMBINED')
            return self.__legal_status

    @property
    def economic_status(self):
        r'''Dictionary of economic status indicators for the chemical.

        Examples
        --------
        >>> pprint(Chemical('benzene').economic_status)
        ["US public: {'Manufactured': 6165232.1, 'Imported': 463146.474, 'Exported': 271908.252}",
         u'1,000,000 - 10,000,000 tonnes per annum',
         u'Intermediate Use Only',
         'OECD HPV Chemicals']
        '''
        if self.__economic_status:
            return self.__economic_status
        else:
            self.__economic_status = economic_status(self.CAS, Method='Combined')
            return self.__economic_status


    @property
    def UNIFAC_groups(self):
        r'''Dictionary of UNIFAC subgroup: count groups for the original
        UNIFAC subgroups, as determined by `DDBST's online service <http://www.ddbst.com/unifacga.html>`_.

        Examples
        --------
        >>> pprint(Chemical('Cumene').UNIFAC_groups)
        {1: 2, 9: 5, 13: 1}
        '''
        if self.__UNIFAC_groups:
            return self.__UNIFAC_groups
        else:
            load_group_assignments_DDBST()
            if self.InChI_Key in DDBST_UNIFAC_assignments:
                self.__UNIFAC_groups = DDBST_UNIFAC_assignments[self.InChI_Key]
                return self.__UNIFAC_groups
            else:
                return None

    @property
    def UNIFAC_Dortmund_groups(self):
        r'''Dictionary of Dortmund UNIFAC subgroup: count groups for the
        Dortmund UNIFAC subgroups, as determined by `DDBST's online service <http://www.ddbst.com/unifacga.html>`_.

        Examples
        --------
        >>> pprint(Chemical('Cumene').UNIFAC_Dortmund_groups)
        {1: 2, 9: 5, 13: 1}
        '''
        if self.__UNIFAC_Dortmund_groups:
            return self.__UNIFAC_Dortmund_groups
        else:
            load_group_assignments_DDBST()
            if self.InChI_Key in DDBST_MODIFIED_UNIFAC_assignments:
                self.__UNIFAC_Dortmund_groups = DDBST_MODIFIED_UNIFAC_assignments[self.InChI_Key]
                return self.__UNIFAC_Dortmund_groups
            else:
                return None

    @property
    def PSRK_groups(self):
        r'''Dictionary of PSRK subgroup: count groups for the PSRK subgroups,
        as determined by `DDBST's online service <http://www.ddbst.com/unifacga.html>`_.

        Examples
        --------
        >>> pprint(Chemical('Cumene').PSRK_groups)
        {1: 2, 9: 5, 13: 1}
        '''
        if self.__PSRK_groups:
            return self.__PSRK_groups
        else:
            load_group_assignments_DDBST()
            if self.InChI_Key in DDBST_PSRK_assignments:
                self.__PSRK_groups = DDBST_PSRK_assignments[self.InChI_Key]
                return self.__PSRK_groups
            else:
                return None

    @property
    def UNIFAC_R(self):
        r'''UNIFAC `R` (normalized Van der Waals volume), dimensionless.
        Used in the UNIFAC model.

        Examples
        --------
        >>> Chemical('benzene').UNIFAC_R
        3.1878
        '''
        if self.UNIFAC_groups:
            return UNIFAC_RQ(self.UNIFAC_groups)[0]
        return None

    @property
    def UNIFAC_Q(self):
        r'''UNIFAC `Q` (normalized Van der Waals area), dimensionless.
        Used in the UNIFAC model.

        Examples
        --------
        >>> Chemical('decane').UNIFAC_Q
        6.016
        '''
        if self.UNIFAC_groups:
            return UNIFAC_RQ(self.UNIFAC_groups)[1]
        return None

    @property
    def Van_der_Waals_volume(self):
        r'''Unnormalized Van der Waals volume, in units of [m^3/mol].

        Examples
        --------
        >>> Chemical('hexane').Van_der_Waals_volume
        6.8261966e-05
        '''
        if self.UNIFAC_R:
            return Van_der_Waals_volume(self.UNIFAC_R)
        return None

    @property
    def Van_der_Waals_area(self):
        r'''Unnormalized Van der Waals area, in units of [m^2/mol].

        Examples
        --------
        >>> Chemical('hexane').Van_der_Waals_area
        964000.0
        '''
        if self.UNIFAC_Q:
            return Van_der_Waals_area(self.UNIFAC_Q)
        return None

    ### One phase properties - calculate lazily
    @property
    def Psat(self):
        r'''Vapor pressure of the chemical at its current temperature, in units
        of [Pa]. For calculation of this property at other temperatures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface :obj:`thermo.vapor_pressure.VaporPressure`;
        each Chemical instance creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).Psat
        10533.614271198725
        >>> Chemical('water').VaporPressure.T_dependent_property(320)
        10533.614271198725
        >>> Chemical('water').VaporPressure.all_methods
        set(['VDI_PPDS', 'BOILING_CRITICAL', 'WAGNER_MCGARRY', 'AMBROSE_WALTON', 'COOLPROP', 'LEE_KESLER_PSAT', 'EOS', 'ANTOINE_POLING', 'SANJARI', 'DIPPR_PERRY_8E', 'Edalat'])
        '''
        return self.VaporPressure.evaluate(self.T)

    @property
    def Hvap(self):
        r'''Enthalpy of vaporization of the chemical at its current temperature,
        in units of [J/mol]. For calculation of this property at other
        temperatures, or specifying manually the method used to calculate it,
        and more - see the object oriented interface
        :obj:`thermo.phase_change.EnthalpyVaporization`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).Hvapm
        43048.23612280223
        >>> Chemical('water').EnthalpyVaporization.T_dependent_property(320)
        43048.23612280223
        >>> Chemical('water').EnthalpyVaporization.all_methods
        set(['VDI_PPDS', 'MORGAN_KOBAYASHI', 'VETERE', 'VELASCO', 'LIU', 'COOLPROP', 'CRC_HVAP_298', 'CLAPEYRON', 'SIVARAMAN_MAGEE_KOBAYASHI', 'ALIBAKHSHI', 'DIPPR_PERRY_8E', 'RIEDEL', 'CHEN', 'PITZER', 'CRC_HVAP_TB'])
        '''
        return self.EnthalpyVaporization.evaluate(self.T)


    @property
    def Cps(self):
        r'''Solid-phase heat capacity of the chemical at its current temperature,
        in units of [J/mol/K]. For calculation of this property at other
        temperatures, or specifying manually the method used to calculate it,
        and more - see the object oriented interface
        :obj:`thermo.heat_capacity.HeatCapacitySolid`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('palladium').Cpsm
        24.930765664000003
        >>> Chemical('palladium').HeatCapacitySolid.T_dependent_property(320)
        25.098979200000002
        >>> Chemical('palladium').HeatCapacitySolid.all_methods
        set(["Perry's Table 2-151", 'CRC Standard Thermodynamic Properties of Chemical Substances', 'Lastovka, Fulem, Becerra and Shaw (2008)'])
        '''
        return self.HeatCapacitySolid.evaluate(self.T)

    @property
    def Cpl(self):
        r'''Liquid-phase heat capacity of the chemical at its current temperature,
        in units of [J/mol/K]. For calculation of this property at other
        temperatures, or specifying manually the method used to calculate it,
        and more - see the object oriented interface
        :obj:`thermo.heat_capacity.HeatCapacityLiquid`; each Chemical instance
        creates one to actually perform the calculations.

        Notes
        -----
        Some methods give heat capacity along the saturation line, some at
        1 atm but only up to the normal boiling point, and some give heat
        capacity at 1 atm up to the normal boiling point and then along the
        saturation line. Real-liquid heat capacity is pressure dependent, but
        this interface is not.

        Examples
        --------
        >>> Chemical('water').Cplm
        75.31462591538556
        >>> Chemical('water').HeatCapacityLiquid.T_dependent_property(320)
        75.2591744360631
        >>> Chemical('water').HeatCapacityLiquid.T_dependent_property_integral(300, 320)
        1505.0619005000553
        '''
        return self.HeatCapacityLiquid.evaluate(self.T)

    @property
    def Cpg(self):
        r'''Gas-phase ideal gas heat capacity of the chemical at its current
        temperature, in units of [J/mol/K]. For calculation of this property at
        other temperatures, or specifying manually the method used to calculate
        it, and more - see the object oriented interface
        :obj:`thermo.heat_capacity.HeatCapacityGas`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water').Cpgm
        33.583577868850675
        >>> Chemical('water').HeatCapacityGas.T_dependent_property(320)
        33.67865044005934
        >>> Chemical('water').HeatCapacityGas.T_dependent_property_integral(300, 320)
        672.6480417835064
        '''
        return self.HeatCapacityGas.evaluate(self.T)

    @property
    def Cvg(self):
        r'''Gas-phase ideal-gas contant-volume heat capacity of the chemical at
        its current temperature, in units of [J/mol/K]. Subtracts R from
        the ideal-gas heat capacity; does not include pressure-compensation
        from an equation of state.
        '''
        return self.HeatCapacityGas.evaluate(self.T) - _R

    @property
    def isentropic_exponent(self):
        r'''Gas-phase ideal-gas isentropic exponent of the chemical at its
        current temperature, [dimensionless]. Does not include
        pressure-compensation from an equation of state.

        Examples
        --------
        >>> Chemical('hydrogen').isentropic_exponent
        1.405237786321222
        '''
        return isentropic_exponent(self.Cpg, self.Cvg)

    @property
    def Vs(self):
        r'''Solid-phase molar volume of the chemical at its current
        temperature, in units of [mol/m^3]. For calculation of this property at
        other temperatures, or specifying manually the method used to calculate
        it, and more - see the object oriented interface
        :obj:`thermo.volume.VolumeSolid`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('iron').Vs
        7.09593392630242e-06
        '''
        return self.VolumeSolid.evaluate(self.T, self.P)

    @property
    def Vl(self):
        r'''Liquid-phase molar volume of the chemical at its current
        temperature and pressure, in units of [mol/m^3]. For calculation of this
        property at other temperatures or pressures, or specifying manually the
        method used to calculate it, and more - see the object oriented interface
        :obj:`thermo.volume.VolumeLiquid`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('cyclobutane', T=225).Vl
        7.42395423425395e-05
        '''
        return self.VolumeLiquid.evaluate(self.T, self.P)

    @property
    def Vg(self):
        r'''Gas-phase molar volume of the chemical at its current
        temperature and pressure, in units of [mol/m^3]. For calculation of this
        property at other temperatures or pressures, or specifying manually the
        method used to calculate it, and more - see the object oriented interface
        :obj:`thermo.volume.VolumeGas`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        Estimate the molar volume of the core of the sun, at 15 million K and
        26.5 PetaPascals, assuming pure helium (actually 68% helium):

        >>> Chemical('helium', T=15E6, P=26.5E15).Vg
        4.805464238181197e-07
        '''
        return self.VolumeGas.evaluate(self.T, self.P)

    @property
    def rhos(self):
        r'''Solid-phase mass density of the chemical at its current temperature,
        in units of [kg/m^3]. For calculation of this property at
        other temperatures, or specifying manually the method used
        to calculate it, and more - see the object oriented interface
        :obj:`thermo.volume.VolumeSolid`; each Chemical instance
        creates one to actually perform the calculations. Note that that
        interface provides output in molar units.

        Examples
        --------
        >>> Chemical('iron').rhos
        7869.999999999994
        '''
        return self.Vms * self.MW / 1000.

    @property
    def rhol(self):
        r'''Liquid-phase mass density of the chemical at its current
        temperature and pressure, in units of [kg/m^3]. For calculation of this
        property at other temperatures and pressures, or specifying manually
        the method used to calculate it, and more - see the object oriented
        interface :obj:`thermo.volume.VolumeLiquid`; each Chemical instance
        creates one to actually perform the calculations. Note that that
        interface provides output in molar units.

        Examples
        --------
        >>> Chemical('o-xylene', T=297).rhol
        876.9946785618097
        '''
        return self.Vml * self.MW / 1000.

    @property
    def rhog(self):
        r'''Gas-phase mass density of the chemical at its current temperature
        and pressure, in units of [kg/m^3]. For calculation of this property at
        other temperatures or pressures, or specifying manually the method used
        to calculate it, and more - see the object oriented interface
        :obj:`thermo.volume.VolumeGas`; each Chemical instance
        creates one to actually perform the calculations. Note that that
        interface provides output in molar units.

        Examples
        --------
        Estimate the density of the core of the sun, at 15 million K and
        26.5 PetaPascals, assuming pure helium (actually 68% helium):

        >>> Chemical('helium', T=15E6, P=26.5E15).rhog
        8329.27226509739

        Compared to a result on
        `Wikipedia <https://en.wikipedia.org/wiki/Solar_core>`_ of 150000
        kg/m^3, the fundamental equation of state performs poorly.

        >>> He = Chemical('helium', T=15E6, P=26.5E15)
        >>> He.VolumeGas.set_user_methods_P(['IDEAL']); He.rhog
        850477.8065477367

        The ideal-gas law performs somewhat better, but vastly overshoots
        the density prediction.
        '''
        return self.Vmg * self.MW /1000.

    @property
    def Zs(self):
        r'''Compressibility factor of the chemical in the solid phase at the
        current temperature and pressure, [dimensionless].

        Utilizes the object oriented interface and
        :obj:`thermo.volume.VolumeSolid` to perform the actual calculation of
        molar volume.

        Examples
        --------
        >>> Chemical('palladium').Z
        0.00036248477437931853
        '''
        return Z(self.T, self.P, self.Vs)

    @property
    def Zl(self):
        r'''Compressibility factor of the chemical in the liquid phase at the
        current temperature and pressure, [dimensionless].

        Utilizes the object oriented interface and
        :obj:`thermo.volume.VolumeLiquid` to perform the actual calculation of
        molar volume.

        Examples
        --------
        >>> Chemical('water').Zl
        0.0007385375470263454
        '''
        return Z(self.T, self.P, self.Vl)

    @property
    def Zg(self):
        r'''Compressibility factor of the chemical in the gas phase at the
        current temperature and pressure, [dimensionless].

        Utilizes the object oriented interface and
        :obj:`thermo.volume.VolumeGas` to perform the actual calculation of
        molar volume.

        Examples
        --------
        >>> Chemical('sulfur hexafluoride', T=700, P=1E9).Zg
        11.140084184207813
        '''
        return Z(self.T, self.P, self.Vg)
        
    @property
    def Bvirial(self):
        r'''Second virial coefficient of the gas phase of the chemical at its
        current temperature and pressure, in units of [mol/m^3].

        This property uses the object-oriented interface
        :obj:`thermo.volume.VolumeGas`, converting its result with
        :obj:`thermo.utils.B_from_Z`.

        Examples
        --------
        >>> Chemical('water').Bvirial
        -0.0009596286322838357
        '''
        return B_from_Z(self.Zg, self.T, self.P)

    @property
    def isobaric_expansion_l(self):
        r'''Isobaric (constant-pressure) expansion of the liquid phase of the
        chemical at its current temperature and pressure, in units of [1/K].

        .. math::
            \beta = \frac{1}{V}\left(\frac{\partial V}{\partial T} \right)_P

        Utilizes the temperature-derivative method of
        :obj:`thermo.volume.VolumeLiquid` to perform the actual calculation.
        The derivatives are all numerical.

        Examples
        --------
        >>> Chemical('dodecane', T=400).isobaric_expansion_l
        0.0011617555762469477
        '''
        dV_dT = self.VolumeLiquid.differentiate_T(self.T, self.P)
        return isobaric_expansion(V=self.Vl, dV_dT=dV_dT)

    @property
    def isobaric_expansion_g(self):
        r'''Isobaric (constant-pressure) expansion of the gas phase of the
        chemical at its current temperature and pressure, in units of [1/K].

        .. math::
            \beta = \frac{1}{V}\left(\frac{\partial V}{\partial T} \right)_P

        Utilizes the temperature-derivative method of
        :obj:`thermo.VolumeGas` to perform the actual calculation.
        The derivatives are all numerical.

        Examples
        --------
        >>> Chemical('Hexachlorobenzene', T=900).isobaric_expansion_g
        0.001151869741981048
        '''
        dV_dT = self.VolumeGas.differentiate_T(self.T, self.P)
        return isobaric_expansion(V=self.Vg, dV_dT=dV_dT)

    @property
    def mul(self):
        r'''Viscosity of the chemical in the liquid phase at its current
        temperature and pressure, in units of [Pa*s].

        For calculation of this property at other temperatures and pressures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface
        :obj:`thermo.viscosity.ViscosityLiquid`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).mul
        0.0005767262693751547
        '''
        return self.ViscosityLiquid.evaluate(self.T, self.P)

    @property
    def mug(self):
        r'''Viscosity of the chemical in the gas phase at its current
        temperature and pressure, in units of [Pa*s].

        For calculation of this property at other temperatures and pressures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface
        :obj:`thermo.viscosity.ViscosityGas`; each Chemical instance
        creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320, P=100).mug
        1.0431450856297212e-05
        '''
        return self.ViscosityGas.evaluate(self.T, self.P)

    @property
    def kl(self):
        r'''Thermal conductivity of the chemical in the liquid phase at its
        current temperature and pressure, in units of [W/m/K].

        For calculation of this property at other temperatures and pressures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface
        :obj:`thermo.thermal_conductivity.ThermalConductivityLiquid`; each
        Chemical instance creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).kl
        0.6369957248212118
        '''
        return self.ThermalConductivityLiquid.evaluate(self.T, self.P)

    @property
    def kg(self):
        r'''Thermal conductivity of the chemical in the gas phase at its
        current temperature and pressure, in units of [W/m/K].

        For calculation of this property at other temperatures and pressures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface
        :obj:`thermo.thermal_conductivity.ThermalConductivityGas`; each
        Chemical instance creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).kg
        0.021273128263091207
        '''
        return self.ThermalConductivityGas.evaluate(self.T, self.P)

    @property
    def sigma(self):
        r'''Surface tension of the chemical at its current temperature, in
        units of [N/m].

        For calculation of this property at other temperatures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface :obj:`thermo.interface.SurfaceTension`;
        each Chemical instance creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('water', T=320).sigma
        0.06855002575793023
        >>> Chemical('water', T=320).SurfaceTension.solve_prop(0.05)
        416.8307110842183
        '''
        return self.SurfaceTension.evaluate(self.T)

    @property
    def permittivity(self):
        r'''Relative permittivity (dielectric constant) of the chemical at its 
        current temperature, [dimensionless].

        For calculation of this property at other temperatures,
        or specifying manually the method used to calculate it, and more - see
        the object oriented interface :obj:`thermo.permittivity.Permittivity`;
        each Chemical instance creates one to actually perform the calculations.

        Examples
        --------
        >>> Chemical('toluene', T=250).permittivity
        2.49775625
        '''
        return self.Permittivity.evaluate(self.T)
    
    @property
    def absolute_permittivity(self):
        r'''Absolute permittivity of the chemical at its current temperature,
        in units of [farad/meter]. Those units are equivalent to 
        ampere^2*second^4/kg/m^3.

        Examples
        --------
        >>> Chemical('water', T=293.15).absolute_permittivity
        7.096684821859018e-10
        '''
        return self.permittivity*physical_constants['electric constant'][0]

    @property
    def JTl(self):
        r'''Joule Thomson coefficient of the chemical in the liquid phase at
        its current temperature and pressure, in units of [K/Pa].

        .. math::
            \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{C_p}
            \left[T \left(\frac{\partial V}{\partial T}\right)_P - V\right]
            = \frac{V}{C_p}\left(\beta T-1\right)

        Utilizes the temperature-derivative method of
        :obj:`thermo.volume.VolumeLiquid` and the temperature-dependent heat
        capacity method :obj:`thermo.heat_capacity.HeatCapacityLiquid` to
        obtain the properties required for the actual calculation.

        Examples
        --------
        >>> Chemical('dodecane', T=400).JTl
        -3.0827160465192742e-07
        '''
        return Joule_Thomson(T=self.T, V=self.Vl, Cp=self.Cplm, beta=self.isobaric_expansion_l)

    @property
    def JTg(self):
        r'''Joule Thomson coefficient of the chemical in the gas phase at
        its current temperature and pressure, in units of [K/Pa].

        .. math::
            \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{C_p}
            \left[T \left(\frac{\partial V}{\partial T}\right)_P - V\right]
            = \frac{V}{C_p}\left(\beta T-1\right)

        Utilizes the temperature-derivative method of
        :obj:`thermo.volume.VolumeGas` and the temperature-dependent heat
        capacity method :obj:`thermo.heat_capacity.HeatCapacityGas` to
        obtain the properties required for the actual calculation.

        Examples
        --------
        >>> Chemical('dodecane', T=400, P=1000).JTg
        5.4089897835384913e-05
        '''
        return Joule_Thomson(T=self.T, V=self.Vg, Cp=self.Cpgm, beta=self.isobaric_expansion_g)

    @property
    def nul(self):
        r'''Kinematic viscosity of the liquid phase of the chemical at its
        current temperature and pressure, in units of [m^2/s].

        .. math::
            \nu = \frac{\mu}{\rho}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.volume.VolumeLiquid`,
        :obj:`thermo.viscosity.ViscosityLiquid`  to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('methane', T=110).nul
        2.858088468937331e-07
        '''
        return nu_mu_converter(mu=self.mul, rho=self.rhol)

    @property
    def nug(self):
        r'''Kinematic viscosity of the gas phase of the chemical at its
        current temperature and pressure, in units of [m^2/s].

        .. math::
            \nu = \frac{\mu}{\rho}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.volume.VolumeGas`,
        :obj:`thermo.viscosity.ViscosityGas`  to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('methane', T=115).nug
        2.5056924327995865e-06
        '''
        return nu_mu_converter(mu=self.mug, rho=self.rhog)

    @property
    def alphal(self):
        r'''Thermal diffusivity of the liquid phase of the chemical at its
        current temperature and pressure, in units of [m^2/s].

        .. math::
            \alpha = \frac{k}{\rho Cp}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.volume.VolumeLiquid`,
        :obj:`thermo.thermal_conductivity.ThermalConductivityLiquid`,
        and :obj:`thermo.heat_capacity.HeatCapacityLiquid` to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('nitrogen', T=70).alphal
        9.444949636299626e-08
        '''
        return thermal_diffusivity(k=self.kl, rho=self.rhol, Cp=self.Cpl)

    @property
    def alphag(self):
        r'''Thermal diffusivity of the gas phase of the chemical at its
        current temperature and pressure, in units of [m^2/s].

        .. math::
            \alpha = \frac{k}{\rho Cp}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.volume.VolumeGas`,
        :obj:`thermo.thermal_conductivity.ThermalConductivityGas`,
        and :obj:`thermo.heat_capacity.HeatCapacityGas` to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('ammonia').alphag
        1.6931865425158556e-05
        '''
        return thermal_diffusivity(k=self.kg, rho=self.rhog, Cp=self.Cpg)

    @property
    def Prl(self):
        r'''Prandtl number of the liquid phase of the chemical at its
        current temperature and pressure, [dimensionless].

        .. math::
            Pr = \frac{C_p \mu}{k}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.viscosity.ViscosityLiquid`,
        :obj:`thermo.thermal_conductivity.ThermalConductivityLiquid`,
        and :obj:`thermo.heat_capacity.HeatCapacityLiquid` to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('nitrogen', T=70).Prl
        2.7828214501488886
        '''
        return Prandtl(Cp=self.Cpl, mu=self.mul, k=self.kl)

    @property
    def Prg(self):
        r'''Prandtl number of the gas phase of the chemical at its
        current temperature and pressure, [dimensionless].

        .. math::
            Pr = \frac{C_p \mu}{k}

        Utilizes the temperature and pressure dependent object oriented
        interfaces :obj:`thermo.viscosity.ViscosityGas`,
        :obj:`thermo.thermal_conductivity.ThermalConductivityGas`,
        and :obj:`thermo.heat_capacity.HeatCapacityGas` to calculate the
        actual properties.

        Examples
        --------
        >>> Chemical('NH3').Prg
        0.847263731933008
        '''
        return Prandtl(Cp=self.Cpg, mu=self.mug, k=self.kg)

    @property
    def Parachor(self):
        r'''Parachor of the chemical at its current temperature and pressure, in units of [N^0.25*m^2.75/mol].

        .. math::
            P = \frac{\sigma^{0.25} MW}{\rho_L - \rho_V}

        Calculated based on surface tension, density of the liquid and gas
        phase, and molecular weight. For uses of this property, see
        :obj:`thermo.utils.Parachor`.

        Examples
        --------
        >>> Chemical('octane').Parachor
        6.291693072841486e-05
        '''
        return Parachor(sigma=self.sigma, MW=self.MW, rhol=self.rhol, rhog=self.rhog)

    ### Single-phase properties
    @phase_property
    def Cp(self):
        r'''Molar heat capacity of the chemical at its current phase and
        temperature, in units of [J/mol/K].

        Utilizes the object oriented interfaces
        :obj:`thermo.heat_capacity.HeatCapacitySolid`,
        :obj:`thermo.heat_capacity.HeatCapacityLiquid`,
        and :obj:`thermo.heat_capacity.HeatCapacityGas` to perform the
        actual calculation of each property.

        Examples
        --------
        >>> Chemical('cubane').Cp
        137.05489206785944
        >>> Chemical('ethylbenzene', T=550, P=3E6).Cp
        294.18449553310046
        '''

    @phase_property
    def V(self):
        r'''Molar volume of the chemical at its current phase and
        temperature and pressure, in units of [m^3/mol].

        Utilizes the object oriented interfaces
        :obj:`thermo.volume.VolumeSolid`,
        :obj:`thermo.volume.VolumeLiquid`,
        and :obj:`thermo.volume.VolumeGas` to perform the
        actual calculation of each property.

        Examples
        --------
        >>> Chemical('ethylbenzene', T=550, P=3E6).V
        0.00017758024401627633
        '''

    @property
    def Z(self):
        r'''Compressibility factor of the chemical at its current phase and
        temperature and pressure, [dimensionless].

        Examples
        --------
        >>> Chemical('MTBE', T=900, P=1E-2).Z
        0.9999999999079768
        '''
        return Z(self.T, self.P, self.V)
    
    @phase_property
    def isobaric_expansion(self):
        r'''Isobaric (constant-pressure) expansion of the chemical at its
        current phase and temperature, in units of [1/K].

        .. math::
            \beta = \frac{1}{V}\left(\frac{\partial V}{\partial T} \right)_P

        Examples
        --------
        Radical change  in value just above and below the critical temperature
        of water:

        >>> Chemical('water', T=647.1, P=22048320.0).isobaric_expansion
        0.34074205839222449

        >>> Chemical('water', T=647.2, P=22048320.0).isobaric_expansion
        0.18143324022215077
        '''

    @phase_property
    def JT(self):
        r'''Joule Thomson coefficient of the chemical at its
        current phase and temperature, in units of [K/Pa].

        .. math::
            \mu_{JT} = \left(\frac{\partial T}{\partial P}\right)_H = \frac{1}{C_p}
            \left[T \left(\frac{\partial V}{\partial T}\right)_P - V\right]
            = \frac{V}{C_p}\left(\beta T-1\right)

        Examples
        --------
        >>> Chemical('water').JT
        -2.2150394958666407e-07
        '''

    @phase_property
    def mu(self):
        r'''Viscosity of the chemical at its current phase, temperature, and
        pressure in units of [Pa*s].

        Utilizes the object oriented interfaces
        :obj:`thermo.viscosity.ViscosityLiquid` and
        :obj:`thermo.viscosity.ViscosityGas` to perform the
        actual calculation of each property.

        Examples
        --------
        >>> Chemical('ethanol', T=300).mu
        0.001044526538460911
        >>> Chemical('ethanol', T=400).mu
        1.1853097849748217e-05
        '''

    @phase_property
    def k(self):
        r'''Thermal conductivity of the chemical at its current phase,
        temperature, and pressure in units of [W/m/K].

        Utilizes the object oriented interfaces
        :obj:`thermo.thermal_conductivity.ThermalConductivityLiquid` and
        :obj:`thermo.thermal_conductivity.ThermalConductivityGas` to perform
        the actual calculation of each property.

        Examples
        --------
        >>> Chemical('ethanol', T=300).kl
        0.16313594741877802
        >>> Chemical('ethanol', T=400).kg
        0.026019924109310026
        '''

    @phase_property
    def nu(self):
        r'''Kinematic viscosity of the the chemical at its current temperature,
        pressure, and phase in units of [m^2/s].

        .. math::
            \nu = \frac{\mu}{\rho}

        Examples
        --------
        >>> Chemical('argon').nu
        1.3846930410865003e-05
        '''

    @phase_property
    def alpha(self):
        r'''Thermal diffusivity of the chemical at its current temperature,
        pressure, and phase in units of [m^2/s].

        .. math::
            \alpha = \frac{k}{\rho Cp}

        Examples
        --------
        >>> Chemical('furfural').alpha
        8.696537158635412e-08
        '''

    @phase_property
    def Pr(self):
        r'''Prandtl number of the chemical at its current temperature,
        pressure, and phase; [dimensionless].

        .. math::
            Pr = \frac{C_p \mu}{k}

        Examples
        --------
        >>> Chemical('acetone').Pr
        4.183039103542709
        '''

    @property
    def Poynting(self):
        r'''Poynting correction factor [dimensionless] for use in phase 
        equilibria methods based on activity coefficients or other reference 
        states. Performs the shortcut calculation assuming molar volume is 
        independent of pressure.

        .. math::
            \text{Poy} =  \exp\left[\frac{V_l (P-P^{sat})}{RT}\right]

        The full calculation normally returns values very close to the
        approximate ones. This property is defined in terms of
        pure components only.

        Examples
        --------
        >>> Chemical('pentane', T=300, P=1E7).Poynting
        1.5743051250679803

        Notes
        -----
        The full equation shown below can be used as follows:

        .. math::
            \text{Poy} = \exp\left[\frac{\int_{P_i^{sat}}^P V_i^l dP}{RT}\right]

        >>> from scipy.integrate import quad
        >>> c = Chemical('pentane', T=300, P=1E7)
        >>> exp(quad(lambda P : c.VolumeLiquid(c.T, P), c.Psat, c.P)[0]/R/c.T)
        1.5821826990975127
        '''
        return exp(self.Vl*(self.P-self.Psat)/_R/self.T)

    def Tsat(self, P):
        return self.VaporPressure.solve(P)

    ### Convenience Dimensionless numbers
    def Reynolds(self, V=None, D=None):
        return Reynolds(V=V, D=D, rho=self.rho, mu=self.mu)

    def Capillary(self, V=None):
        return Capillary(V=V, mu=self.mu, sigma=self.sigma)

    def Weber(self, V=None, D=None):
        return Weber(V=V, L=D, rho=self.rho, sigma=self.sigma)

    def Bond(self, L=None):
        return Bond(rhol=self.rhol, rhog=self.rhog, sigma=self.sigma, L=L)

    def Jakob(self, Tw=None):
        return Jakob(Cp=self.Cp, Hvap=self.Hvap, Te=Tw-self.T)

    def Grashof(self, Tw=None, L=None):
        return Grashof(L=L, beta=self.isobaric_expansion, T1=Tw, T2=self.T,
                       rho=self.rho, mu=self.mu)

    def Peclet_heat(self, V=None, D=None):
        return Peclet_heat(V=V, L=D, rho=self.rho, Cp=self.Cp, k=self.k)

