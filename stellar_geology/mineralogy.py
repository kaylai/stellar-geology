"""CIPW normative mineralogy calculations.

Provides functions to compute normative mineral assemblages (olivine,
clinopyroxene, orthopyroxene, garnet) from silicate compositions using the
Thompson (1982) reaction-space approach as implemented in Putirka and
Rarick (2019).
"""

import numpy as np
import pandas as pd

from . import constants as const
from . import conversions as conv

# ---------------------------------------------------------------------------
# Cation-balance matrix: rows = minerals, columns = oxide components
# Cells DR6:DU9 in Putirka and Rarick (2019) supplemental spreadsheet 2.
# ---------------------------------------------------------------------------
MINERAL_TRANSFORMATIONS = pd.DataFrame(
    [[1,   0,   2,   0  ],
     [1.8, 0.1, 1.4, 0.6],
     [1.8, 0.1, 1.8, 0.2],
     [3,   1,   2.7, 0.3]],
    index=["olivine", "clinopyroxene", "orthopyroxene", "garnet"],
    columns=["SiO2", "Al2O3", "FmO", "CaO"]
)

MINERAL_TRANSFORMATIONS_INV = pd.DataFrame(
    np.linalg.inv(MINERAL_TRANSFORMATIONS.values),
    index=MINERAL_TRANSFORMATIONS.columns,
    columns=MINERAL_TRANSFORMATIONS.index
)

def calculate_mineralogy(silicate_composition: dict[str, float], units: str = 'wtpt_oxides',
                         pressure=1, model="PutirkaRarick19") -> dict[str, float]:
    """
    Calculates the CIPW normative mineralogy for olivine (ol), clinopyroxene (cpx), orthopyroxene
    (opx), and garnet (gar). Uses equations of Thompson (1982) "Reaction space: An algebraic and
    geometric approach" as implemented in Putirka and Rarick (2019) supplemental spreadsheet 2.

    Parameters
    ----------
    silicate_composition : dict[str, float]
        Dictionary of chemical components describing a bulk silicate composition. Must include, at
        a minimum, keys for SiO2, Al2O3, FeO, MgO, and CaO (any units).
    units : str
        Units of the silicate_composition dict. Defaults to 'wtpt_oxides'. Any valid unit string is
        accepted (e.g., 'wtpt_elements', 'molfrac_oxides'). Compositions are converted to wt% oxides
        internally.
    pressure : float
        Only required if using model of Till2012. Pressure in GPa. Default is 1.
    model : str
        Which modal mineralogy calculation to use:
            - PutirkaRarick19: Putirka and Rarick (2019) [Default]
            - Till2012_plag: Till et al. 2012 JGR. Plagioclase lherzolite model
            - Till2012_sp: Till et al. 2012 JGR. Spinel lherzolite model

    Returns
    -------
    dict[str, float]
        Normative mineralogy with keys "olivine", "clinopyroxene", "orthopyroxene", "garnet".
    
    # TODO: Implement mineralogy class like the VESIcal Calculate() object and subclasses that
    then get passed onto either P&R or Till2012 in its own file. This will make it easier to 
    expand the library to add more mineralogy calculation models.
    """
    import warnings as w
    
    if units not in conv.VALID_UNITS:
        raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

    # Convert inputs to canonical wtpt_oxides
    if silicate_composition is not None and units != 'wtpt_oxides':
        silicate_composition = conv.convert_to_wtpt_oxides(silicate_composition, units)
    else:
        silicate_composition = conv.normalize_composition(silicate_composition, 'wtpt_oxides')

    # Fill missing CIPW oxides with 0 and warn
    if model=="PutirkaRarick19":
        cipw_oxides = ["SiO2", "Al2O3", "FeO", "MgO", "CaO"]
    elif model=="Till2012_plag" or model=="Till2012_sp":
        cipw_oxides = ["TiO2", "FeO", "MgO", "CaO", "Na2O", "K2O"]
    else:
        raise ValueError(f"{model} not a recognized model.")
        
    missing = [ox for ox in cipw_oxides if ox not in silicate_composition]
    if missing:
        w.warn(f"{missing} missing from silicate composition and will be "
            "set to 0 for mineralogy calculation.", category=UserWarning)
        for ox in missing:
            silicate_composition[ox] = 0.0

    # Route to individual models
    if model=="PutirkaRarick19":
        mol_prop_ox_cipw = _calculate_mol_prop_ox_cipw(silicate_composition)
        mol_frac_cipw = _calculate_mol_frac_cipw(mol_prop_ox_cipw)
    
        # excel "MMULT" on inverted matrix and mol_frac_cipw's
        mol_frac_cipw_series = pd.Series(mol_frac_cipw)
        result = mol_frac_cipw_series.reindex(MINERAL_TRANSFORMATIONS_INV.index) @ MINERAL_TRANSFORMATIONS_INV
        
        return dict(result)
    
    if model=="Till2012_plag":
        return _calculate_mineralogy_Till(composition=silicate_composition,
                                pressure=pressure,
                                regression="plagioclase")
    if model=="Till2012_sp":
        return _calculate_mineralogy_Till(composition=silicate_composition,
                                pressure=pressure,
                                regression="spinel")

# ---------------------------------------------------------------------------
# PUTIRKA
# ---------------------------------------------------------------------------  

def _calculate_mol_prop_ox_cipw(silicate_composition: dict[str, float]) -> dict[str, float]:
    xtal_oxides = ["SiO2", "Al2O3", "FeO", "MgO", "CaO"]
    mol_prop_ox_cipw = {ox: (silicate_composition.get(ox, 0.0) * const.CationNum[ox]) /
                        const.oxideMass[ox] for ox in xtal_oxides}
    return mol_prop_ox_cipw

def _calculate_mol_frac_cipw(mol_prop_ox_cipw: dict[str, float]) -> dict[str, float]:
    mpoc = mol_prop_ox_cipw
    # direct copy for Si, Al, Ca
    mol_frac_cipw = {k: v/sum(list(mpoc.values())) for k, v in mpoc.items() if k not in ["FeO", "MgO"]}
    # FmO = FeO + MgO
    mol_frac_cipw["FmO"] = (mpoc["FeO"] + mpoc["MgO"])/sum(list(mpoc.values()))

    return mol_frac_cipw

# ---------------------------------------------------------------------------
# TILL
# --------------------------------------------------------------------------- 

def _calculate_mineralogy_Till(composition: dict[str, float], pressure=1, regression="plagioclase") -> dict[str, float]:
    """
    Uses model equation of Till et al. 2012 for plagioclase lherzolite.
    Model form: mode_mineral = a + b*P + c*(1-Mg#) + d*NaK# + e*TiO2 + f*K2O
        - where mode_mineral is modal abundance of mineral (either olivine, cpx, plag, or qtz),
        P is pressure in GPa, Mg# is molar Mg/(Mg+Fe), NaK# is wt% (Na2O+K2O)/(Na2O+K2O+CaO), TiO2 
        and K2O are in wt%, and a, b, c, d, e, and f are constants.
    
    | Mineral | a | b | c | d | e | f |
    |---------|-----|-----|-----|-----|-----|-----|
    | **Plag Lherz** | | | | | | |
    | Olivine | 0.132 | 0.069 | 0.173 | -0.261 | -0.009 | 0.010 |
    | Clinopyroxene | 0.238 | -0.051 | -0.065 | -0.097 | 0.004 | -0.013 |
    | Plagioclase | 0.408 | 0.146 | -0.136 | 0.630 | 0.014 | -0.028 |
    | Quartz | 0.222 | -0.164 | -0.029 | -0.272 | -0.009 | 0.031 |
    | **Sp Lherz** | | | | | | |
    | Olivine | 0.123 | 0.085 | 0.203 | -0.268 | -0.019 | 0.012 |
    | Clinopyroxene | 0.180 | 0.005 | 0.075 | -0.264 | 0.007 | 0.012 |
    | Plagioclase | 0.554 | -0.043 | -0.192 | 0.949 | -0.007 | -0.041 |
    | Quartz | 0.143 | -0.046 | -0.085 | -0.416 | 0.020 | 0.018 |
    
    Parameters
    ----------
    composition : dict[str, float]
        Dictionary of oxide composition of a bulk silicate planet (or mantle) in wt%.
    pressure : float
        Pressure of the calculation, in GPa. Default is 1 GPa.
    regression : str
        Which set of model coefficients to use. "plagioclase" uses Plagioclase Lherzolite
        regressed coefficients. "spinel" uses Spinel Lherzolite regressed coefficients. See
        Table 5 of Till et al. (2012) JGR.
    """
    
    mg_number = conv.calculate_mg_number(composition)
    NaK_number = (
        (composition['Na2O'] + composition['K2O']) / 
        (composition['Na2O'] + composition['K2O'] + composition['CaO'])
    )
    
    if regression == "plagioclase":
        # TABLE 5 OF PUBLICATION
        # Version 4 regressions in Christy's spreadsheets
        # test currently for "Plag Lherz All Data" tab, unknown version of regression
        mineral_constants = {
            "olivine" : {
                'a' : 0.132, 'b' : 0.069, 'c' : 0.173, 'd' : -0.261, 'e' : -0.009, 'f' : 0.010,},
            "clinopyroxene" : {
                'a' : 0.238, 'b' : -0.051, 'c' : -0.065, 'd' : -0.097, 'e' : 0.004, 'f' : -0.013,},
            "plagioclase" : {
                'a' : 0.408, 'b' : 0.146, 'c' : -0.136, 'd' : 0.630, 'e' : 0.014, 'f' : -0.028,},
            "quartz" : {
                'a' : 0.222, 'b' : -0.164, 'c' : -0.029, 'd' : -0.272, 'e' : -0.009, 'f' : 0.031,},
        }
        
    elif regression == "spinel":
        # TABLE 5 OF PUBLICATION
        mineral_constants = {
            "olivine" : {
                'a' : 0.123, 'b' : 0.085, 'c' : 0.203, 'd' : -0.268, 'e' : -0.019, 'f' : 0.012,},
            "clinopyroxene" : {
                'a' : 0.180, 'b' : 0.005, 'c' : 0.075, 'd' : -0.264, 'e' : 0.007, 'f' : 0.012,},
            "plagioclase" : {
                'a' : 0.554, 'b' : -0.043, 'c' : -0.192, 'd' : 0.949, 'e' : -0.007, 'f' : -0.041,},
            "quartz" : {
                'a' : 0.143, 'b' : -0.046, 'c' : -0.085, 'd' : -0.416, 'e' : 0.020, 'f' : 0.018,},
        }
    else:
        raise ValueError(f"{regression} is not a recognized model regression. Must be 'plagioclase'"
                          "or 'spinel'.")
    
    mineralogy_dict = {}
    for mineral, v in mineral_constants.items():
        calc_mode = (v['a'] + v['b']*pressure + v['c']*(1-mg_number) +
                     v['d']*NaK_number + v['e']*composition['TiO2'] +
                     v['f']*composition['K2O'])
        mineralogy_dict[mineral] = calc_mode
    
    
    # TESTING ONES IN REGRESSION CALCS THAT MATCH "ALL DATA" TAB IN SPREADSHEET
    # Version 6 in Christy's spreadsheets - not completed updating - check vals
    # mineral_constants = {
    #     "olivine" : {
    #         'a' : 0.130, 'b' : 0.074, 'c' : 0.160, 'd' : -0.266, 'e' : -0.007, 'f' : 0.014,},
    #     "clinopyroxene" : {
    #         'a' : 0.238, 'b' : -0.051, 'c' : -0.065, 'd' : -0.097, 'e' : 0.004, 'f' : -0.013,},
    #     "plagioclase" : {
    #         'a' : 0.408, 'b' : 0.146, 'c' : -0.136, 'd' : 0.630, 'e' : 0.014, 'f' : -0.028,},
    #     "quartz" : {
    #         'a' : 0.222, 'b' : -0.164, 'c' : -0.029, 'd' : -0.272, 'e' : -0.009, 'f' : 0.031,},
    # }
    
    mineralogy_normed = {
        mineral: mode/sum(mineralogy_dict.values()) for mineral, mode in mineralogy_dict.items()}
    
    return mineralogy_normed
      
# ---------------------------------------------------------------------------
# REVERSE PIPELINE CALCULATIONS (CURRENTLY PUTIRKA ONLY)
# ---------------------------------------------------------------------------       

def _mol_frac_cipw_to_wtpt_oxides(mol_prop_ox_cipw: dict[str, float]) -> dict[str, float]:
    """
    Converts CIPW mol fractions back to wt% oxides. This is the reverse of
    _calculate_mol_prop_ox_cipw followed by normalization.

    Forward:  mol_prop[ox] = wt%[ox] * CationNum[ox] / oxideMass[ox]
    Reverse:  wt_raw[ox]   = mol_frac[ox] * oxideMass[ox] / CationNum[ox]
    Then normalize to sum to 100.

    Parameters
    ----------
    mol_prop_ox_cipw : dict[str, float]
        Dictionary of oxide mol proportions with keys SiO2, Al2O3, FeO, MgO, CaO.
        Values need not sum to 1; they are renormalized internally.

    Returns
    -------
    dict[str, float]
        Composition as wt% oxides normalized to 100.
    """
    wt_raw = {ox: mf * const.oxideMass[ox] / const.CationNum[ox]
              for ox, mf in mol_prop_ox_cipw.items()}
    wt_sum = sum(wt_raw.values())
    return {ox: 100.0 * v / wt_sum for ox, v in wt_raw.items()}

def calculate_composition_from_mineralogy(mineralogy: dict[str, float], mg_number: float) -> dict[str, float]:
    """
    Recovers a bulk silicate planet composition (wt% oxides) from CIPW normative mineralogy.
    This is the partial inverse of calculate_mineralogy().

    Because the forward pipeline compresses FeO + MgO into a single FmO component,
    an assumed Mg# is required to split FmO back into FeO and MgO.

    Only the 5 CIPW oxides (SiO2, Al2O3, FeO, MgO, CaO) are recovered; minor
    oxides (TiO2, Na2O, Cr2O3, NiO) discarded in the forward pipeline cannot be
    reconstructed.

    Parameters
    ----------
    mineralogy : dict[str, float]
        Dictionary with keys "olivine", "clinopyroxene", "orthopyroxene", "garnet"
        and float values (mol fractions, as returned by calculate_mineralogy).
    mg_number : float
        Molar Mg# = Mg/(Mg+Fe), range (0, 1]. Controls how the mafic oxides mineralogy component,
        FmO, is split into FeO and MgO in the reverse pipeline.

    Returns
    -------
    dict[str, float]
        Bulk silicate composition as wt% oxides (SiO2, Al2O3, FeO, MgO, CaO),
        normalized to 100.
    """
    required_phases = ["olivine", "clinopyroxene", "orthopyroxene", "garnet"]
    missing = [p for p in required_phases if p not in mineralogy]
    if missing:
        raise ValueError(f"Missing required mineral phases: {missing}")

    if not (0 < mg_number <= 1):
        raise ValueError(f"mg_number must be in (0, 1], got {mg_number}")

    mineralogy_series = pd.Series(mineralogy).reindex(MINERAL_TRANSFORMATIONS.index)
    mol_frac_cipw = mineralogy_series @ MINERAL_TRANSFORMATIONS

    fmo = mol_frac_cipw["FmO"]
    mol_prop_ox_cipw = {
        "SiO2":  mol_frac_cipw["SiO2"],
        "Al2O3": mol_frac_cipw["Al2O3"],
        "FeO":   fmo * (1.0 - mg_number),
        "MgO":   fmo * mg_number,
        "CaO":   mol_frac_cipw["CaO"],
    }

    return _mol_frac_cipw_to_wtpt_oxides(mol_prop_ox_cipw)

def plot_norm(mineralogy: dict[str, float]) -> dict[str, float]:
    """
    Normalizes mineralogy to only olivine, clinopyroxene, and orthopyroxene for ternary plotting.
    
    Parameters
    ----------
    mineralogy : dict[str, float]
        Dictionary with key:value pairs as mineral name: proportion. Requires at minimum "olivine",
        "clinopyroxene", and "orthopyroxene". All other key:value pairs will be ignored. Put
        whatever you want there, I don't care. But you won't get it returned to you.

    Returns
    -------
    dict[str, float]
        Mineralogy normalized to ol, opx, and cpx only.
    """
    required_phases = ["olivine", "clinopyroxene", "orthopyroxene"]
    missing = [p for p in required_phases if p not in mineralogy]
    if missing:
        raise ValueError(f"Missing required mineral phases: {missing}")

    normed_phases = {k: v for k, v in mineralogy.items() if k in required_phases}
    sum_phases = sum(normed_phases.values())

    return {k: v / sum_phases for k, v in normed_phases.items()}

# class Mineralogy(object):
#     def __init__(self, bulk_silicate_planet, units='wtpt_oxides'):
#         """
#         Returns a Planet() object.

#         Parameters
#         ----------
#         bulk_silicate_planet : dict
#             Bulk silicate planet composition. Units specified by the `units`
#             parameter.
#         units : str
#             Units of the bulk_planet and bulk_silicate_planet dicts. Defaults
#             to 'wtpt_oxides'. Any valid unit string is accepted (e.g.,
#             'wtpt_elements', 'molfrac_oxides'). Compositions are converted
#             to wt% oxides internally on init.

#         Returns
#         -------
#         Mineralogy() object.
#         """
#         if units not in conv.VALID_UNITS:
#             raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")
        
#         # Convert inputs to canonical wtpt_oxides
#         if bulk_silicate_planet is not None and units != 'wtpt_oxides':
#             bulk_silicate_planet = conv.convert_to_wtpt_oxides(bulk_silicate_planet, units)
        
#         self._bulk_silicate_planet = bulk_silicate_planet
#         self._mineralogy = None

#     @property
#     def bulk_silicate_planet(self):
#         if self._bulk_silicate_planet is not None: # unneeded but good catch for future implem.
#             return self._bulk_silicate_planet
    
#     @property
#     def mineralogy(self):
#         if self._mineralogy is not None:
#             pass
    
#     @classmethod
#     def from_planet(cls, planet):
#         return cls(bulk_silicate_planet=planet.mineralogy)
