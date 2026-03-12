from . import constants as const
from . import conversions as conv
import numpy as np
import pandas as pd

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
    

def calculate_mineralogy(silicate_composition, units='wtpt_oxides'):
    """
    Calculates the CIPW normative mineralogy for olivine (ol), clinopyroxene (cpx), orthpyroxene
    (opx), and garnet (gar). Uses equations of Thompson (1982) "Reaction space: An algebraic and
    geometric approach" as implemented in Putirka and Rarick (2019) supplemental spreadsheet 2.
    
    Parameters
    ----------
    silicate_composition : dict
        Dictionary of chemical components describing a bulk silicate composition. Must include, at
        a minimum, keys for SiO2, Al2O3, FeO, MgO, and CaO (any units)
    units : str
        Units of the silicate_composition dict. Defaults to 'wtpt_oxides'. Any valid unit string is
        accepted (e.g., 'wtpt_elements', 'molfrac_oxides'). Compositions are converted to wt% oxides
        internally.
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
    
    # calculate mineral proportions as the product of two arrays
    # cation sums for each mineral: 3 for ol, 2 for px, and 8 for gt
    # cells DR6:DU9 in Putirka and Rarick (2019) supplemental spreadsheet 2
    mineral_transformations = pd.DataFrame(
        [[1,   0,   2,   0  ],
        [1.8, 0.1, 1.4, 0.6],
        [1.8, 0.1, 1.8, 0.2],
        [3,   1,   2.7, 0.3]],
        index=["olivine", "clinopyroxene", "orthopyroxene", "garnet"],
        columns=["SiO2", "Al2O3", "FmO", "CaO"]
    )
    
    # excel "MINVERSE" on cation sums matrix
    mineral_transformations_inv = pd.DataFrame(
        np.linalg.inv(mineral_transformations.values),
        index=mineral_transformations.columns,       # rows are now oxides
        columns=mineral_transformations.index        # columns are now minerals
    )
    
    # excel "MMULT" on inverted matrix and mol_frac_cipw's
    mol_frac_cipw_series = pd.Series(mol_frac_cipw)
    result = mol_frac_cipw_series.reindex(mineral_transformations_inv.index) @ mineral_transformations_inv
    
    mineralogy = dict(result)
    
    return mineralogy
    
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
    
def plot_norm(mineralogy):
    """
    Normalizes mineralogy to only olivine, clinopyroxene, and orthopyroxene for ternary plotting.
    
    Paramaters
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
    phases = ["olivine", "clinopyroxene", "orthopyroxene"]
    not_found = [k for k in phases if k not in mineralogy.keys()]
    if not_found is not None:
        if len(not_found) > 0:
            raise ValueError(f"{not_found} not found in {mineralogy}.")
    normed_phases = {k: v for k, v in mineralogy.items() if k in phases}
    
    sum_phases = sum(normed_phases.values())
    mineralogy_norm = {k: v/sum_phases for k, v in normed_phases.items()}
    
    return mineralogy_norm
