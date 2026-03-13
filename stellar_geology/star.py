"""The Star module defines the :class:`Star` class for representing stellar
compositions and converting them through the dex-to-oxide pipeline.
"""

from . import conversions as conv
from . import constants as const
import warnings as w

__all__ = ['Star']


class Star(object):
    # TODO: add `units` param to __init__ so users can pass compositions in
    # any supported unit, not just dex. Internally convert to dex on init.
    def __init__(self, stellar_dex: dict[str, float] | None = None, name: str | None = None, mass: float | None = None) -> None:
        """
        Parameters
        ----------
        stellar_dex : dict[str, float], optional
            Stellar composition in dex notation, as elements.
        name : str, optional
            Arbitrary name for your star as a string. Can be anything. Either
            aa, either bb, even zombocom. At zombocom you can do anything.
        mass : float, optional
            Star mass in some units that I don't know because this isn't
            implemented yet. So, put whatever float you want here. It won't
            make any difference.
        """
        self._stellar_dex = stellar_dex
        self._name = name
        self._mass = mass

        # calculated attributes
        # typically not needed by user, used for benchmarking and debugging
        self._ax: dict[str, float] | None = None
        self._atoms_ref_solar: dict[str, float] | None = None
        self._total_wt_atoms: dict[str, float] | None = None
        self._wtpt_elements: dict[str, float] | None = None
        self._wtpt_oxides: dict[str, float] | None = None

        if stellar_dex is not None:
            unrecognized_keys = []
            for k in stellar_dex.keys():
                if k not in const.elements_to_oxides.keys() and k not in const.oxides_to_elements.keys():
                    unrecognized_keys.append(k)
            if len(unrecognized_keys) > 0:
                w.warn(f"{unrecognized_keys} were not recognized as compositional parameters and "
                       "will be ignored in calculations.", category=UserWarning)


    @property
    def stellar_dex(self) -> dict[str, float] | None:
        """dict[str, float] or None : Stellar composition in dex notation."""
        return self._stellar_dex

    @property
    def name(self) -> str | None:
        """str or None : Star name."""
        return self._name

    @property
    def mass(self) -> float | None:
        """float or None : Star mass (not yet used in calculations)."""
        return self._mass

    @property
    def ax(self) -> dict[str, float] | None:
        """dict[str, float] or None : Elemental ratio relative to solar (10^dex)."""
        if self._ax is None and self._stellar_dex is not None:
            self._ax = conv.calculate_ax_from_dex(self._stellar_dex)
        return self._ax

    @property
    def atoms_ref_solar(self) -> dict[str, float] | None:
        """dict[str, float] or None : Number of atoms referenced to solar abundances."""
        if self._atoms_ref_solar is None and self.ax is not None:
            self._atoms_ref_solar = conv.calculate_atoms_ref_solar_from_ax(self.ax)
        return self._atoms_ref_solar

    @property
    def total_wt_atoms(self) -> dict[str, float] | None:
        """dict[str, float] or None : Total weight of atoms (element wt scaled by atomic mass)."""
        if self._total_wt_atoms is None and self.atoms_ref_solar is not None:
            self._total_wt_atoms = conv.calculate_total_wt_atoms_from_atoms_ref_solar(self.atoms_ref_solar)
        return self._total_wt_atoms

    @property
    def wtpt_elements(self) -> dict[str, float] | None:
        """dict[str, float] or None : Composition as wt% elements (includes volatiles C, O, S)."""
        if self._wtpt_elements is None and self.total_wt_atoms is not None:
            self._wtpt_elements = conv.calculate_wtpt_elements_from_total_wt_atoms(self.total_wt_atoms)
        return self._wtpt_elements

    @property
    def wtpt_oxides(self) -> dict[str, float] | None:
        """dict[str, float] or None : Composition as wt% oxides (volatile-free)."""
        if self._wtpt_oxides is None and self.wtpt_elements is not None:
            self._wtpt_oxides = conv.calculate_wtpt_oxides_from_wtpt_elements(self.wtpt_elements)
        return self._wtpt_oxides

    def get_composition(self, units: str = 'wtpt_oxides', normalization: str | None = None) -> dict[str, float] | None:
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
        dict[str, float] or None
            Composition in the requested units, or None if no stellar_dex is set.

        Notes
        -----
        Star element outputs (wtpt_elements, wtfrac_elements) include volatile
        elements (C, O, S) from the stellar pipeline. These are not available
        in oxide-based conversions.
        """
        if self._stellar_dex is None:
            w.warn("Cannot compute composition: no stellar_dex is set.",
                   category=UserWarning)
            return None

        if units not in conv.VALID_UNITS:
            raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

        # Special case: element wt outputs include volatile elements (C, O, S)
        # from the stellar pipeline that aren't present in oxide-based conversions
        if units in ('wtpt_elements', 'wtfrac_elements'):
            wt_elements = self.wtpt_elements
            if wt_elements is None:
                w.warn("Cannot compute element composition from stellar_dex.",
                       category=UserWarning)
                return None
            if units == 'wtpt_elements':
                result = dict(wt_elements)
            else:
                result = {k: v / 100.0 for k, v in wt_elements.items()}
        else:
            wt_oxides = self.wtpt_oxides
            if wt_oxides is None:
                w.warn("Cannot compute oxide composition from stellar_dex.",
                       category=UserWarning)
                return None
            result = conv.convert_composition(wt_oxides, units)

        if normalization is not None and normalization != 'none':
            if units == 'molfrac_singleO':
                w.warn("Normalization is not supported for units='molfrac_singleO' "
                       "and will be ignored.", category=UserWarning)
            else:
                result = conv.normalize_composition(result, normalization, units)

        return result
