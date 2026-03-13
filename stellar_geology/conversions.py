"""
Chemistry conversion functions imported by planet and star
"""

from . import constants as const

def calculate_bulk_planet_from_dex(stellar_dex: dict[str, float]) -> dict[str, float]:
    """
    Runs all intermediate calculations to go directly from star composition in
    dex notation to bulk planet in wt% oxides. This is the most common use case,
    but allows the other intermediate methods to be exposed for benchmarking and
    debugging.

    Parameters
    ----------
    stellar_dex: dict[str, float]
        Star composition in dex notation.

    Returns
    -------
    dict[str, float]
        Bulk planet composition in wt% oxides.
    """
    ax = calculate_ax_from_dex(stellar_dex)
    atomsRefSolar = calculate_atomsRefSolar_from_ax(ax)
    totalWtAtoms = calculate_totalWtAtoms_from_atomsRefSolar(atomsRefSolar)
    wtptElements = calcualte_wtptElements_from_totalWtAtoms(totalWtAtoms)
    wtptOxides = calculate_wtptOxides_from_wtptElements(wtptElements)
    
    return wtptOxides

#--- INTERMEDIATE CALCULATIONS BETWEEN DEX NOTATION AND BULK PLANET OXIDES ---#
def calculate_ax_from_dex(stellar_dex: dict[str, float]) -> dict[str, float]:
    """
    Convert from dex system notation to elemental ratio relative to solar

    Args:
        composition (dict[str, float])
    """
    ax = {}
    dex_elems = list(stellar_dex.keys())
    
    for el in dex_elems:
        if stellar_dex[el] != 0:
            ax[el] = 10**stellar_dex[el]
        else:
            ax[el] = 0
    
    return ax
    
def calculate_atomsRefSolar_from_ax(ax: dict[str, float]) -> dict[str, float]:
    atomsRefSolar = {}
    ax_elems = list(ax.keys())
    
    for el in ax_elems:
        atomsRefSolar[el] = ax[el] * 10**const.A_El[el]
    
    return atomsRefSolar
    
def calculate_totalWtAtoms_from_atomsRefSolar(atomsRefSolar: dict[str, float]) -> dict[str, float]:
    totalWtAtoms = {}
    atomsRefSolar_elems = list(atomsRefSolar.keys())
    
    for el in atomsRefSolar_elems:
        totalWtAtoms[el] = atomsRefSolar[el] * const.cationMass[el]
    
    return totalWtAtoms

def calcualte_wtptElements_from_totalWtAtoms(totalWtAtoms: dict[str, float]) -> dict[str, float]:
    totalWtAtoms_sum = sum(totalWtAtoms.values())
    wtptElements = {}
    totalWtAtoms_elems = list(totalWtAtoms.keys())
    
    for el in totalWtAtoms_elems:
        wtptElements[el] = 100 * totalWtAtoms[el]/totalWtAtoms_sum
    
    return wtptElements

def calculate_wtptOxides_from_wtptElements(wtptElements: dict[str, float]) -> dict[str, float]:
    wtptOxides = {}
    volatile_free_elems = list(const.elements_to_oxides.keys())
    
    for el in volatile_free_elems:
        if el not in wtptElements:
            continue
        ox = const.elements_to_oxides[el]
        conversion_factor = const.oxideMass[ox]/(const.cationMass[el]*const.CationNum[ox])
        wtptOxides[ox] = wtptElements[el]*conversion_factor
    
    wtptOxides_sum = sum(wtptOxides.values())
    wtptOxides = {k: 100*v/wtptOxides_sum for k, v in wtptOxides.items()}

    return wtptOxides


#--- COMPOSABLE UNIT CONVERSION SYSTEM ---#

# Valid unit strings for convert_composition
VALID_UNITS: list[str] = [
    'wtpt_oxides', 'wtpt_elements',
    'wtfrac_oxides', 'wtfrac_elements',
    'molfrac_oxides', 'molfrac_elements', 'molfrac_singleO',
    'molpt_oxides', 'molpt_elements',
]


def _wt_to_mol_oxides(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """Compute raw numerators: wtpt / oxideMass for each oxide.

    Dividing by the sum of the result yields mole fractions.
    """
    return {oxide: wtpt / const.oxideMass[oxide]
            for oxide, wtpt in wtpt_oxides.items()}


def _wt_to_mol_elements(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """Compute raw numerators: CationNum * wtpt / oxideMass for each element.

    Dividing by the sum of the result yields mole fractions of cations.
    """
    raw = {}
    for oxide, wtpt in wtpt_oxides.items():
        element = const.oxides_to_elements[oxide]
        raw[element] = const.CationNum[oxide] * wtpt / const.oxideMass[oxide]
    return raw


def _wt_oxides_to_wt_elements(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """Compute raw numerators: wtpt * cationMass * CationNum / oxideMass.

    Dividing by the sum and multiplying by 100 yields wt% elements.
    """
    raw = {}
    for oxide, wtpt in wtpt_oxides.items():
        element = const.oxides_to_elements[oxide]
        raw[element] = wtpt * const.cationMass[element] * const.CationNum[oxide] / const.oxideMass[oxide]
    return raw


def _wt_to_mol_singleO(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """Convert wt% oxides to moles of cations per single oxygen atom.

    Normalized per oxygen atom, not by sum of cations.
    """
    cation_moles = {}
    total_O = 0.0

    for oxide, wtpt in wtpt_oxides.items():
        element = const.oxides_to_elements[oxide]
        cation_moles[element] = const.CationNum[oxide] * wtpt / const.oxideMass[oxide]
        total_O += const.OxygenNum[oxide] * wtpt / const.oxideMass[oxide]

    if total_O == 0:
        return dict(cation_moles)
    return {k: v / total_O for k, v in cation_moles.items()}


def convert_composition(wtpt_oxides: dict[str, float], units: str) -> dict[str, float]:
    """Convert a composition in wt% oxides to any supported unit system.

    This is the main dispatcher for all unit conversions. It takes a canonical
    wt% oxides dict and converts it to the requested units by:
    1. Computing raw numerators via the appropriate helper
    2. Dividing by the sum to complete the conversion
    3. Scaling to the target (100 for percent, 1.0 for fraction)

    Parameters
    ----------
    wtpt_oxides : dict[str, float]
        Composition in wt% oxides (oxide keys, values summing to ~100).
    units : str
        Target unit string. One of: 'wtpt_oxides', 'wtpt_elements',
        'wtfrac_oxides', 'wtfrac_elements', 'molfrac_oxides',
        'molfrac_elements', 'molfrac_singleO', 'molpt_oxides',
        'molpt_elements'.

    Returns
    -------
    dict[str, float]
        Composition in the requested units.

    Raises
    ------
    ValueError
        If units is not a recognized unit string.
    """
    if units not in VALID_UNITS:
        raise ValueError(f"units must be one of {VALID_UNITS}, got '{units}'.")

    # Special case: molfrac_singleO is normalized per oxygen atom, not by sum
    if units == 'molfrac_singleO':
        return _wt_to_mol_singleO(wtpt_oxides)

    # Determine the raw numerators and target scale
    raw: dict[str, float]
    if units in ('wtpt_oxides', 'wtfrac_oxides'):
        raw = dict(wtpt_oxides)
    elif units in ('wtpt_elements', 'wtfrac_elements'):
        raw = _wt_oxides_to_wt_elements(wtpt_oxides)
    elif units in ('molfrac_oxides', 'molpt_oxides'):
        raw = _wt_to_mol_oxides(wtpt_oxides)
    elif units in ('molfrac_elements', 'molpt_elements'):
        raw = _wt_to_mol_elements(wtpt_oxides)
    else:
        raise ValueError(f"Unhandled units: '{units}'")

    # Divide by sum and scale to target
    raw_sum = sum(raw.values())
    if raw_sum == 0:
        return dict(raw)

    if 'pt' in units:
        target = 100.0
    else:
        target = 1.0

    return {k: target * v / raw_sum for k, v in raw.items()}


def _wt_elements_to_wt_oxides(wt_elements: dict[str, float]) -> dict[str, float]:
    """Convert element weight values to oxide weight values.

    Reverse of _wt_oxides_to_wt_elements. Multiply each element's weight by
    oxideMass / (cationMass * CationNum), then normalize to sum to 100.
    """
    raw = {}
    for element, wt in wt_elements.items():
        oxide = const.elements_to_oxides[element]
        raw[oxide] = wt * const.oxideMass[oxide] / (const.cationMass[element] * const.CationNum[oxide])

    raw_sum = sum(raw.values())
    if raw_sum == 0:
        return dict(raw)
    return {k: 100.0 * v / raw_sum for k, v in raw.items()}


def convert_to_wtpt_oxides(composition: dict[str, float], from_units: str) -> dict[str, float]:
    """Convert a composition from any supported unit system back to wt% oxides.

    This is the inverse of convert_composition(). Takes a dict in any supported
    unit and returns wt% oxides (oxide keys, values summing to 100).

    Parameters
    ----------
    composition : dict[str, float]
        Composition in the units specified by from_units.
    from_units : str
        Unit string describing the input. One of the strings in VALID_UNITS.

    Returns
    -------
    dict[str, float]
        Composition in wt% oxides, normalized to sum to 100.

    Raises
    ------
    ValueError
        If from_units is not a recognized unit string.
    """
    if from_units not in VALID_UNITS:
        raise ValueError(f"from_units must be one of {VALID_UNITS}, got '{from_units}'.")

    if from_units == 'wtpt_oxides':
        return dict(composition)

    if from_units == 'wtfrac_oxides':
        return {k: v * 100.0 for k, v in composition.items()}

    if from_units == 'wtpt_elements':
        return _wt_elements_to_wt_oxides(composition)

    if from_units == 'wtfrac_elements':
        scaled = {k: v * 100.0 for k, v in composition.items()}
        return _wt_elements_to_wt_oxides(scaled)

    if from_units == 'molfrac_oxides':
        return mol_oxides_to_wtpt_oxides(composition)

    if from_units == 'molfrac_elements':
        return mol_cations_to_wtpt_oxides(composition)

    if from_units == 'molpt_oxides':
        frac = {k: v / 100.0 for k, v in composition.items()}
        return mol_oxides_to_wtpt_oxides(frac)

    if from_units == 'molpt_elements':
        frac = {k: v / 100.0 for k, v in composition.items()}
        return mol_cations_to_wtpt_oxides(frac)

    if from_units == 'molfrac_singleO':
        # Normalize cation ratios to sum to 1.0, then treat as mol_cations
        total = sum(composition.values())
        if total == 0:
            return {const.elements_to_oxides[k]: 0.0 for k in composition}
        frac = {k: v / total for k, v in composition.items()}
        return mol_cations_to_wtpt_oxides(frac)

    # Should be unreachable — all VALID_UNITS are handled above
    raise ValueError(f"Unhandled from_units: '{from_units}'")


#--- LEGACY UNIT CONVERSION FUNCTIONS (kept for reverse conversions) ---#
def wtpt_oxides_to_mol_oxides(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """
    Convert from wt% oxides to mol fraction oxides.

    Parameters
    ----------
    wtpt_oxides: dict[str, float]
        Composition in wt% oxides.

    Returns
    -------
    dict[str, float]
        Composition in mol fraction oxides, normalized to sum to 1.0.
    """
    mol_oxides = {}
    for oxide, wtpt in wtpt_oxides.items():
        mol_oxides[oxide] = wtpt / const.oxideMass[oxide]

    mol_sum = sum(mol_oxides.values())
    if mol_sum == 0:
        return dict(mol_oxides)
    mol_oxides = {k: v / mol_sum for k, v in mol_oxides.items()}

    return mol_oxides


def wtpt_oxides_to_mol_cations(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """
    Convert from wt% oxides to mol fraction cations.

    Parameters
    ----------
    wtpt_oxides: dict[str, float]
        Composition in wt% oxides.

    Returns
    -------
    dict[str, float]
        Composition in mol fraction cations (element keys), normalized to sum
        to 1.0.
    """
    mol_cations = {}
    for oxide, wtpt in wtpt_oxides.items():
        element = const.oxides_to_elements[oxide]
        mol_cations[element] = const.CationNum[oxide] * wtpt / const.oxideMass[oxide]

    mol_sum = sum(mol_cations.values())
    if mol_sum == 0:
        return dict(mol_cations)
    mol_cations = {k: v / mol_sum for k, v in mol_cations.items()}

    return mol_cations


def wtpt_oxides_to_mol_singleO(wtpt_oxides: dict[str, float]) -> dict[str, float]:
    """
    Convert from wt% oxides to moles of cations per single oxygen atom.

    Parameters
    ----------
    wtpt_oxides: dict[str, float]
        Composition in wt% oxides.

    Returns
    -------
    dict[str, float]
        Composition as cation moles normalized to one oxygen atom (element
        keys). Not normalized to sum to 1.0.
    """
    cation_moles = {}
    total_O = 0.0

    for oxide, wtpt in wtpt_oxides.items():
        element = const.oxides_to_elements[oxide]
        cation_moles[element] = const.CationNum[oxide] * wtpt / const.oxideMass[oxide]
        total_O += const.OxygenNum[oxide] * wtpt / const.oxideMass[oxide]

    if total_O == 0:
        return dict(cation_moles)
    mol_singleO = {k: v / total_O for k, v in cation_moles.items()}

    return mol_singleO


def mol_oxides_to_wtpt_oxides(mol_oxides: dict[str, float]) -> dict[str, float]:
    """
    Convert from mol fraction oxides to wt% oxides.

    Parameters
    ----------
    mol_oxides: dict[str, float]
        Composition in mol fraction oxides.

    Returns
    -------
    dict[str, float]
        Composition in wt% oxides, normalized to sum to 100.
    """
    wtpt_oxides = {}
    for oxide, mol_frac in mol_oxides.items():
        wtpt_oxides[oxide] = mol_frac * const.oxideMass[oxide]

    wtpt_sum = sum(wtpt_oxides.values())
    if wtpt_sum == 0:
        return dict(wtpt_oxides)
    wtpt_oxides = {k: 100 * v / wtpt_sum for k, v in wtpt_oxides.items()}

    return wtpt_oxides


def mol_cations_to_wtpt_oxides(mol_cations: dict[str, float]) -> dict[str, float]:
    """
    Convert from mol fraction cations to wt% oxides.

    Parameters
    ----------
    mol_cations: dict[str, float]
        Composition in mol fraction cations (element keys).

    Returns
    -------
    dict[str, float]
        Composition in wt% oxides, normalized to sum to 100.
    """
    wtpt_oxides = {}
    for element, mol_frac in mol_cations.items():
        oxide = const.elements_to_oxides[element]
        wtpt_oxides[oxide] = mol_frac / const.CationNum[oxide] * const.oxideMass[oxide]

    wtpt_sum = sum(wtpt_oxides.values())
    if wtpt_sum == 0:
        return dict(wtpt_oxides)
    wtpt_oxides = {k: 100 * v / wtpt_sum for k, v in wtpt_oxides.items()}

    return wtpt_oxides


def mol_oxides_to_mol_cations(mol_oxides: dict[str, float]) -> dict[str, float]:
    """
    Convert from mol fraction oxides to mol fraction cations.

    Parameters
    ----------
    mol_oxides: dict[str, float]
        Composition in mol fraction oxides.

    Returns
    -------
    dict[str, float]
        Composition in mol fraction cations (element keys), normalized to sum
        to 1.0.
    """
    mol_cations = {}
    for oxide, mol_frac in mol_oxides.items():
        element = const.oxides_to_elements[oxide]
        mol_cations[element] = mol_frac * const.CationNum[oxide]

    mol_sum = sum(mol_cations.values())
    if mol_sum == 0:
        return dict(mol_cations)
    mol_cations = {k: v / mol_sum for k, v in mol_cations.items()}

    return mol_cations


def mol_cations_to_mol_oxides(mol_cations: dict[str, float]) -> dict[str, float]:
    """
    Convert from mol fraction cations to mol fraction oxides.

    Parameters
    ----------
    mol_cations: dict[str, float]
        Composition in mol fraction cations (element keys).

    Returns
    -------
    dict[str, float]
        Composition in mol fraction oxides, normalized to sum to 1.0.
    """
    mol_oxides = {}
    for element, mol_frac in mol_cations.items():
        oxide = const.elements_to_oxides[element]
        mol_oxides[oxide] = mol_frac / const.CationNum[oxide]

    mol_sum = sum(mol_oxides.values())
    if mol_sum == 0:
        return dict(mol_oxides)
    mol_oxides = {k: v / mol_sum for k, v in mol_oxides.items()}

    return mol_oxides


def normalize_composition(composition: dict[str, float], normalization: str, units: str) -> dict[str, float]:
    """
    Normalize a composition dict according to the specified normalization scheme.

    Parameters
    ----------
    composition: dict[str, float]
        Composition to normalize.
    normalization: str
        One of 'none', 'standard', 'fixedvolatiles', 'additionalvolatiles'.
    units: str
        Any valid unit string from VALID_UNITS (e.g., 'wtpt_oxides',
        'molfrac_oxides', 'molpt_elements'). Used to determine the
        normalization target: 100 for percent units, 1.0 for fraction units.

    Returns
    -------
    dict[str, float]
        Normalized composition.
    """
    valid_normalizations = ['none', 'standard', 'fixedvolatiles', 'additionalvolatiles']

    if normalization not in valid_normalizations:
        raise ValueError(f"normalization must be one of {valid_normalizations}, got '{normalization}'.")
    if units not in VALID_UNITS:
        raise ValueError(f"units must be one of {VALID_UNITS}, got '{units}'.")

    result = dict(composition)

    if normalization == 'none':
        return result

    volatiles = ['H2O', 'CO2']
    target = 100.0 if 'pt' in units else 1.0

    if normalization == 'standard':
        total = sum(result.values())
        if total == 0:
            return result
        result = {k: target * v / total for k, v in result.items()}
        return result

    if normalization == 'fixedvolatiles':
        volatile_sum = sum(result.get(v, 0) for v in volatiles)
        non_volatile_target = target - volatile_sum
        non_volatile_sum = sum(v for k, v in result.items() if k not in volatiles)
        if non_volatile_sum == 0:
            return result
        for k in result:
            if k not in volatiles:
                result[k] = non_volatile_target * result[k] / non_volatile_sum
        return result

    if normalization == 'additionalvolatiles':
        non_volatile_sum = sum(v for k, v in result.items() if k not in volatiles)
        if non_volatile_sum == 0:
            return result
        for k in result:
            if k not in volatiles:
                result[k] = target * result[k] / non_volatile_sum
        return result

    # Should be unreachable — all valid normalizations are handled above
    raise ValueError(f"Unhandled normalization: '{normalization}'")
