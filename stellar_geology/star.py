from . import conversions as conv
from . import constants as const
import warnings as w

class Star(object):
    # TODO: add `units` param to __init__ so users can pass compositions in
    # any supported unit, not just dex. Internally convert to dex on init.
    def __init__(self, stellar_dex=None, name=None, mass=None):
        """
        None required. Would return an empty Star() object with no attributes.

        Parameters
        ----------
        stellar_dex:    dict
            Stellar composition in dex notation, as elements.
        name: str
            Arbitrary name for your planet as a string. Can be anything. Either
            aa, either bb, even zombocom. At zombocom you can do anything.
        mass:   float
            Planet mass in some units that I don't know because this isn't
            implemented yet. So, put whatever float you want here. It won't
            make any difference. 
        """
        self._stellar_dex = stellar_dex
        self._name = name
        self._mass = mass
        
        # calculated attributes
        # typically not needed by user, used for benchmarking and debugging
        self._ax = None
        self._atomsRefSolar = None
        self._totalWtAtoms = None
        self._wtptElements = None
        self._wtptOxides = None
        
        if stellar_dex is not None:
            unrecognized_keys = []
            for k in stellar_dex.keys():
                if k not in const.elements_to_oxides.keys() and k not in const.oxides_to_elements.keys():
                    unrecognized_keys.append(k)
            if len(unrecognized_keys) > 0:
                w.warn(f"{unrecognized_keys} were not recognized as compositional parameters and "
                       "will be ignored in calculations.", category=UserWarning)
    
    
    @property
    def stellar_dex(self):
        if self._stellar_dex is not None:
            return self._stellar_dex
    
    @property
    def name(self):
        if self._name is not None:
            return self._name
    
    @property
    def mass(self):
        if self._mass is not None:
            return self._mass
    
    @property
    def ax(self):
        if self._ax is not None:
            return self._ax
        if self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
            return self._ax 
    
    @property
    def atomsRefSolar(self):
        if self._atomsRefSolar is not None:
            return self._atomsRefSolar
        if self._ax is not None:
            return conv.calculate_atomsRefSolar_from_ax(self._ax)
        if self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
            return conv.calculate_atomsRefSolar_from_ax(self._ax)
    
    @property
    def totalWtAtoms(self):
        if self._totalWtAtoms is not None:
            return self._totalWtAtoms
        if self._atomsRefSolar is not None:
            return conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
        if self._ax is not None:
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            return conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
        if self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            return conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
    
    @property
    def wtptElements(self):
        if self._wtptElements is not None:
            return self._wtptElements
        if self._totalWtAtoms is not None:
            return conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
        if self._atomsRefSolar is not None:
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            return conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
        if self._ax is not None:
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            return conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
        if self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            return conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
    
    @property
    def wtptOxides(self):
        if self._wtptOxides is not None:
            return self._wtptOxides
        if self._wtptElements is not None:
            return conv.calculate_wtptOxides_from_wtptElements(self._wtptElements)
        if self._totalWtAtoms is not None:
            self._wtptElements = conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
            return conv.calculate_wtptOxides_from_wtptElements(self._wtptElements) 
        if self._atomsRefSolar is not None:
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            self._wtptElements = conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
            return conv.calculate_wtptOxides_from_wtptElements(self._wtptElements)
        if self._ax is not None:
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            self._wtptElements = conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
            return conv.calculate_wtptOxides_from_wtptElements(self._wtptElements)
        if self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
            self._atomsRefSolar = conv.calculate_atomsRefSolar_from_ax(self._ax)
            self._totalWtAtoms = conv.calculate_totalWtAtoms_from_atomsRefSolar(self._atomsRefSolar)
            self._wtptElements = conv.calcualte_wtptElements_from_totalWtAtoms(self._totalWtAtoms)
            return conv.calculate_wtptOxides_from_wtptElements(self._wtptElements)

    def get_composition(self, units='wtpt_oxides', normalization=None):
        """
        Return the star's composition in the requested units with optional
        normalization.

        Parameters
        ----------
        units : str
            One of 'wtpt_oxides', 'wtpt_elements', 'wtfrac_oxides',
            'wtfrac_elements', 'molfrac_oxides', 'molfrac_elements',
            'molfrac_singleO', 'molpt_oxides', 'molpt_elements'.
        normalization : str or None
            One of None, 'none', 'standard', 'fixedvolatiles',
            'additionalvolatiles'.

        Returns
        -------
        dict or None
            Composition in the requested units, or None if no stellar_dex is set.

        Notes
        -----
        Star element outputs (wtpt_elements, wtfrac_elements) include volatile
        elements (C, O, S) from the stellar pipeline. These are not available
        in oxide-based conversions.
        """
        if self._stellar_dex is None:
            return None

        if units not in conv.VALID_UNITS:
            raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

        # Special case: element wt outputs include volatile elements (C, O, S)
        # from the stellar pipeline that aren't present in oxide-based conversions
        if units == 'wtpt_elements':
            result = dict(self.wtptElements)
        elif units == 'wtfrac_elements':
            result = {k: v / 100.0 for k, v in self.wtptElements.items()}
        else:
            result = conv.convert_composition(self.wtptOxides, units)

        if normalization is not None and normalization != 'none':
            if units == 'molfrac_singleO':
                w.warn("Normalization is not supported for units='molfrac_singleO' "
                       "and will be ignored.", category=UserWarning)
            else:
                result = conv.normalize_composition(result, normalization, units)

        return result
