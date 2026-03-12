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

def calculate_mineralogy(silicate_composition, units='wtpt_oxides'):
    """
    Calculates the CIPW normative mineralogy for olivine (ol), clinopyroxene (cpx), orthopyroxene
    (opx), and garnet (gar). Uses equations of Thompson (1982) "Reaction space: An algebraic and
    geometric approach" as implemented in Putirka and Rarick (2019) supplemental spreadsheet 2.

    Parameters
    ----------
    silicate_composition : dict
        Dictionary of chemical components describing a bulk silicate composition. Must include, at
        a minimum, keys for SiO2, Al2O3, FeO, MgO, and CaO (any units).
    units : str
        Units of the silicate_composition dict. Defaults to 'wtpt_oxides'. Any valid unit string is
        accepted (e.g., 'wtpt_elements', 'molfrac_oxides'). Compositions are converted to wt% oxides
        internally.

    Returns
    -------
    dict
        Normative mineralogy with keys "olivine", "clinopyroxene", "orthopyroxene", "garnet".
    """
    if units not in conv.VALID_UNITS:
        raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

    # Convert inputs to canonical wtpt_oxides
    if silicate_composition is not None and units != 'wtpt_oxides':
        silicate_composition = conv.convert_to_wtpt_oxides(silicate_composition, units)
    else:
        silicate_composition = conv.normalize_composition(silicate_composition, 'standard', 'wtpt_oxides')
    
    mol_prop_ox_cipw = _calculate_mol_prop_ox_cipw(silicate_composition)
    mol_frac_cipw = _calculate_mol_frac_cipw(mol_prop_ox_cipw)
    
    # excel "MMULT" on inverted matrix and mol_frac_cipw's
    mol_frac_cipw_series = pd.Series(mol_frac_cipw)
    result = mol_frac_cipw_series.reindex(MINERAL_TRANSFORMATIONS_INV.index) @ MINERAL_TRANSFORMATIONS_INV
    
    return dict(result)
    
def _calculate_mol_prop_ox_cipw(silicate_composition):
    xtal_oxides = ["SiO2", "Al2O3", "FeO", "MgO", "CaO"]
    mol_prop_ox_cipw = {ox: (silicate_composition[ox] * const.CationNum[ox]) /
                        const.oxideMass[ox] for ox in xtal_oxides}
    return mol_prop_ox_cipw

def _calculate_mol_frac_cipw(mol_prop_ox_cipw):
    mpoc = mol_prop_ox_cipw
    # direct copy for Si, Al, Ca
    mol_frac_cipw = {k: v/sum(list(mpoc.values())) for k, v in mpoc.items() if k not in ["FeO", "MgO"]}
    # FmO = FeO + MgO
    mol_frac_cipw["FmO"] = (mpoc["FeO"] + mpoc["MgO"])/sum(list(mpoc.values()))

    return mol_frac_cipw

def _mol_frac_cipw_to_wtpt_oxides(mol_prop_ox_cipw):
    """
    Converts CIPW mol fractions back to wt% oxides. This is the reverse of
    _calculate_mol_prop_ox_cipw followed by normalization.

    Forward:  mol_prop[ox] = wt%[ox] * CationNum[ox] / oxideMass[ox]
    Reverse:  wt_raw[ox]   = mol_frac[ox] * oxideMass[ox] / CationNum[ox]
    Then normalize to sum to 100.

    Parameters
    ----------
    mol_prop_ox_cipw : dict
        Dictionary of oxide mol proportions with keys SiO2, Al2O3, FeO, MgO, CaO.
        Values need not sum to 1; they are renormalized internally.

    Returns
    -------
    dict
        Composition as wt% oxides normalized to 100.
    """
    wt_raw = {ox: mf * const.oxideMass[ox] / const.CationNum[ox]
              for ox, mf in mol_prop_ox_cipw.items()}
    wt_sum = sum(wt_raw.values())
    return {ox: 100.0 * v / wt_sum for ox, v in wt_raw.items()}

def calculate_composition_from_mineralogy(mineralogy, mg_number=0.89):
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
    mineralogy : dict
        Dictionary with keys "olivine", "clinopyroxene", "orthopyroxene", "garnet"
        and float values (mol fractions, as returned by calculate_mineralogy).
    mg_number : float
        Molar Mg# = Mg/(Mg+Fe), range (0, 1]. Default 0.89 (approximate Earth
        mantle value). Controls the FeO/MgO split in the recovered composition.

    Returns
    -------
    dict
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

def plot_norm(mineralogy):
    """
    Normalizes mineralogy to only olivine, clinopyroxene, and orthopyroxene for ternary plotting.
    
    Parameters
    ----------
    mineralogy : dict[str:float]
        Dictionary with key:value pairs as mineral name: proportion. Requires at minimum "olivine",
        "clinopyroxene", and "orthopyroxene". All other key:value pairs will be ignored. Put 
        whatever you want there, I don't care. But you won't get it returned to you.
    
    Returns
    -------
    dict
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
