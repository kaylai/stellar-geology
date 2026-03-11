"""
Chemistry conversion functions imported by planet and star
"""

from . import constants as const

def calculate_bulk_planet_from_dex(stellar_dex):
    """
    Runs all intermediate calculations to go directly from star composition in
    dex notation to bulk planet in wt% oxides. This is the most common use case,
    but allows the other intermediate methods to be exposed for benchmarking and
    debugging.
    
    Parameters
    ----------
    stellar_dex: dict
        Star composition in dex notation.
        
    Returns
    -------
    dict
        Bulk planet composition in wt% oxides.
    """
    ax = calculate_ax_from_dex(stellar_dex)
    atomsRefSolar = calculate_atomsRefSolar_from_ax(ax)
    totalWtAtoms = calculate_totalWtAtoms_from_atomsRefSolar(atomsRefSolar)
    wtptElements = calcualte_wtptElements_from_totalWtAtoms(totalWtAtoms)
    wtptOxides = calculate_wtptOxides_from_wtptElements(wtptElements)
    
    return wtptOxides

#--- CALCULATIONS BETWEEN BULK PLANET AND BULK SILICATE PLANET ---#
def calculate_silicate_from_bulk(bulk_planet, alphaFe):
    """
    Calculates the bulk silicate composition given known bulk composition and
    the alphaFe ratio for partitioning bulk Fe between the core and mantle, as
    defined in Putirka and Rarick (2019).
    
    Parameters
    ----------
    bulk_planet:    dict
        Bulk planet composition in wt% oxides. Must contain value for FeO.
    alphaFe:    float
            Ratio of Fe in the bulk silicate planet and bulk planet, defined
            in Putirka and Rarick (2019): alphaFe = FeBSP/FeBP. Will always
            be a positive fraction <1.
    
    Returns
    -------
    dict
        Bulk silicate planet composition in wt% oxides.
    """
    if "FeO" not in list(bulk_planet.keys()):
        raise ValueError("Bulk planet composition must have FeO concentration.")
    
    if alphaFe <= 0 or alphaFe >= 1:
        raise ValueError("alphaFe must be a float where 0 < alphaFe < 1")
    
    Fe_bulk = bulk_planet["FeO"]
    Fe_silicate = Fe_bulk * alphaFe
    
    
    pass

def calculate_bulk_from_silicate(bulk_silicate_planet, alphaFe):
    """
    Calculates the bulk planet composition given known bulk silicate composition
    and the alphaFe ratio for partitioning bulk Fe between the core and mantle,
    as defined in Putirka and Rarick (2019).
    
    Parameters
    ----------
    bulk_silicate_planet:    dict
        Bulk silicate planet composition in wt% oxides.
    alphaFe:    float
        Ratio of Fe in the bulk silicate planet and bulk planet, defined
        in Putirka and Rarick (2019): alphaFe = FeBSP/FeBP. Will always
        be a positive fraction <1.
    
    Returns
    -------
    dict
        Bulk planet composition in wt% oxides.
    """
    return alphaFe * bulk_planet

#--- INTERMEDIATE CALCULATIONS BETWEEN DEX NOTATION AND BULK PLANET OXIDES ---#
def calculate_ax_from_dex(stellar_dex):
    """
    Convert from dex system notation to elemental ratio relative to solar

    Args:
        composition (dict)
    """
    ax = {}
    dex_elems = list(stellar_dex.keys())
    
    for el in dex_elems:
        if stellar_dex[el] != 0:
            ax[el] = 10**stellar_dex[el]
        else:
            ax[el] = 0
    
    return ax
    
def calculate_atomsRefSolar_from_ax(ax):
    atomsRefSolar = {}
    ax_elems = list(ax.keys())
    
    for el in ax_elems:
        atomsRefSolar[el] = ax[el] * 10**const.A_El[el]
    
    return atomsRefSolar
    
def calculate_totalWtAtoms_from_atomsRefSolar(atomsRefSolar):
    totalWtAtoms = {}
    atomsRefSolar_elems = list(atomsRefSolar.keys())
    
    for el in atomsRefSolar_elems:
        totalWtAtoms[el] = atomsRefSolar[el] * const.cationMass[el]
    
    return totalWtAtoms

def calcualte_wtptElements_from_totalWtAtoms(totalWtAtoms):
    totalWtAtoms_sum = sum(totalWtAtoms.values())
    wtptElements = {}
    totalWtAtoms_elems = list(totalWtAtoms.keys())
    
    for el in totalWtAtoms_elems:
        wtptElements[el] = 100 * totalWtAtoms[el]/totalWtAtoms_sum
    
    return wtptElements

def calculate_wtptOxides_from_wtptElements(wtptElements):
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
