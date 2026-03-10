from . import constants as const

class Star(object):
    def __init__(self, composition_dex):
        """Describes the chemistry of a star and its rocky planets (bulk, volatile-free)

        Args:
            composition_dex (dict): Stellar composition in dex notation.
        """
        self._composition_dex = composition_dex
    
    def get_composition(self, units="dex"):
        """Return dict of compositional values

        Args:
            units (str): what units to return. Intermediate units are not so useful
            for science but accessible to benchmark calculations against the Putirka
            and Rarick (2019) supplementary data files. Can be one of:
                - "dex" (default)
                - "AX"
                - "atomsRefSolar"
                - "totalWtAtoms"
                - "wtptElements"
                - "wtptOxides"
        """
        comp_dex = self._composition_dex
        
        if units == "dex":
            return comp_dex
        
        comp_AX = self.dex_to_AX(comp_dex)
        if units == "AX":
            return comp_AX
        comp_atomsRefSolar = self.AX_to_atomsRefSolar(comp_AX)
        if units == "atomsRefSolar":
            return comp_atomsRefSolar
        comp_totalWtAtoms = self.atomsRefSolar_to_totalWtAtoms(comp_atomsRefSolar)
        if units == "totalWtAtoms":
            return comp_totalWtAtoms
        comp_wtptElements = self.totalWtAtoms_to_wtptElements(comp_totalWtAtoms)
        if units == "wtptElements":
            return comp_wtptElements
        comp_wtptOxides = self.wtptElements_to_wtptOxides(comp_wtptElements)
        if units == "wtptOxides":
            return comp_wtptOxides
        else:
            raise ValueError("Units must be one of 'dex', 'AX', 'atomsRefSolar', "
                 "'totalWtAtoms', 'wtptElements', or 'wtptOxides'.")

    def dex_to_AX(self, composition):
        """
        Convert from dex system notation to elemental ratio relative to solar

        Args:
            composition (dict)
        """
        AX = {}
        dex_elems = list(composition.keys())
        
        for el in dex_elems:
            if composition[el] != 0:
                AX[el] = 10**composition[el]
            else:
                AX[el] = 0
        
        return AX

    def AX_to_atomsRefSolar(self, composition):
        atomsRefSolar = {}
        AX_elems = list(composition.keys())
        
        for el in AX_elems:
            atomsRefSolar[el] = composition[el] * 10**const.A_El[el]
        
        return atomsRefSolar

    def atomsRefSolar_to_totalWtAtoms(self, composition):
        totalWtAtoms = {}
        atomsRefSolar_elems = list(composition.keys())
        
        for el in atomsRefSolar_elems:
            totalWtAtoms[el] = composition[el] * const.cationMass[el]
        
        return totalWtAtoms

    def totalWtAtoms_to_wtptElements(self, composition):
        totalWtAtoms_sum = sum(composition.values())
        wtptElements = {}
        totalWtAtoms_elems = list(composition.keys())
        
        for el in totalWtAtoms_elems:
            wtptElements[el] = 100 * composition[el]/totalWtAtoms_sum
        
        return wtptElements

    def wtptElements_to_wtptOxides(self, composition):
        wtptOxides = {}
        volatile_free_elems = list(const.elements_to_oxides.keys())
        
        for el in volatile_free_elems:
            if el not in composition:
                continue
            ox = const.elements_to_oxides[el]
            conversion_factor = const.oxideMass[ox]/(const.cationMass[el]*const.CationNum[ox])
            wtptOxides[ox] = composition[el]*conversion_factor
        
        wtptOxides_sum = sum(wtptOxides.values())
        wtptOxides = {k: 100*v/wtptOxides_sum for k, v in wtptOxides.items()}
        
        return wtptOxides
    