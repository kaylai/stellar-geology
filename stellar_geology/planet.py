"""The Planet module defines the :class:`Planet` class for representing
planetary bulk and silicate compositions derived from stellar data or
direct geochemical inputs.
"""

from __future__ import annotations

from typing import Any, TYPE_CHECKING

from . import conversions as conv
from . import constants as const

if TYPE_CHECKING:
    from .star import Star

__all__ = ['Planet']


class Planet(object):
    def __init__(self, bulk_planet: dict[str, float] | None = None,
                 bulk_silicate_planet: dict[str, float] | None = None,
                 stellar_dex: dict[str, float] | None = None,
                 alphas: dict[str, float] | None = None,
                 name: str | None = None, mass: float | None = None,
                 mineralogy: Any = None, units: str = 'wtpt_oxides') -> None:
        """
        Returns a Planet() object.

        Parameters
        ----------
        Caution: do not pass multiple conflicting compositional parameters or
        it will raise an Error. Just pass one, and the others will be auto-
        matically computed for you. Like magic.

        bulk_planet : dict[str, float]
            Bulk planet composition. Units specified by the `units` parameter.
        bulk_silicate_planet : dict[str, float]
            Bulk silicate planet composition. Units specified by the `units`
            parameter.
        stellar_dex : dict[str, float]
            Star composition in dex notation.
        alphas : dict[str, float]
            Ratio of element in the bulk silicate planet and bulk planet, defined
            in Putirka and Rarick (2019): e.g., alphas = FeBSP/FeBP. Will always
            be a positive fraction <1. Used for defining which elements partition
            into a metallic core. Commonly, Fe, Si, and Ni. Fe is required when
            passing this argument: {'Fe': 0.49}.
        name : str
            Arbitrary name for your planet as a string. Can be anything. Either
            aa, either bb, even zombocom. At zombocom you can do anything.
        mass : float
            Planet mass in some units that I don't know because this isn't
            implemented yet. So, put whatever float you want here. It won't
            make any difference.
        mineralogy : Mineralogy() object
            Not yet implemented to do anything if this is input. Must be
            generated with calculate_mineralogy(). In future will be stored
            as Planet attr.
        units : str
            Units of the bulk_planet and bulk_silicate_planet dicts. Defaults
            to 'wtpt_oxides'. Any valid unit string is accepted (e.g.,
            'wtpt_elements', 'molfrac_oxides'). Compositions are converted
            to wt% oxides internally on init.

        Returns
        -------
        Planet() object.
        """
        if units not in conv.VALID_UNITS:
            raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

        if bulk_planet is not None:
            bulk_planet = const.filter_compositional_keys(bulk_planet, 'bulk_planet')
        if bulk_silicate_planet is not None:
            bulk_silicate_planet = const.filter_compositional_keys(
                bulk_silicate_planet, 'bulk_silicate_planet')
        if stellar_dex is not None:
            stellar_dex = const.filter_compositional_keys(stellar_dex, 'stellar_dex')

        if bulk_planet is not None and stellar_dex is not None:
            raise ValueError("Can not pass both bulk_planet and stellar_dex.")
        
        if bulk_planet is not None and bulk_silicate_planet is not None and alphas is not None:
            raise ValueError("Cannot pass all bulk_planet, bulk_silicate_planet"
                             ", and alphas as values may be contradictory.")
        
        if bulk_silicate_planet is not None and alphas is not None and stellar_dex is not None:
            raise ValueError("Cannot pass all bulk_silicate_planet, alphas, "
                             "and stellar_dex, as values may be contradictory.")
        
        # Convert inputs to canonical wtpt_oxides
        if bulk_planet is not None and units != 'wtpt_oxides':
            bulk_planet = conv.convert_to_wtpt_oxides(bulk_planet, units)
        if bulk_silicate_planet is not None and units != 'wtpt_oxides':
            bulk_silicate_planet = conv.convert_to_wtpt_oxides(
                bulk_silicate_planet, units)

        self._bulk_planet = bulk_planet
        self._bulk_silicate_planet = bulk_silicate_planet
        self._stellar_dex = stellar_dex
        self._alphas = alphas
        self._name = name
        self._mass = mass

    @property
    def bulk_planet(self) -> dict[str, float] | None:
        """dict[str, float] or None : Bulk planet composition in wt% oxides.

        Auto-calculated from ``stellar_dex`` or
        ``bulk_silicate_planet`` + ``alphas`` if not provided directly.
        """
        if self._bulk_planet is not None:
            return self._bulk_planet
        if self._stellar_dex is not None:
            return conv.calculate_bulk_planet_from_dex(self._stellar_dex)
        if self._bulk_silicate_planet is not None and self._alphas is not None:
            return self._calculate_bulk_from_silicate(
                bulk_silicate_planet=self._bulk_silicate_planet, alphas=self._alphas)
        return None

    @property
    def bulk_silicate_planet(self) -> dict[str, float] | None:
        """dict[str, float] or None : Bulk silicate planet composition in wt% oxides.

        Auto-calculated from ``bulk_planet`` + ``alphas`` if not provided
        directly.
        """
        if self._bulk_silicate_planet is not None:
            return self._bulk_silicate_planet
        bulk = self.bulk_planet
        if bulk is not None and self._alphas is not None:
            return self._calculate_silicate_from_bulk(
                bulk_planet=bulk, alphas=self._alphas)
        return None

    @property
    def stellar_dex(self) -> dict[str, float] | None:
        """dict[str, float] or None : Stellar composition in dex notation.

        Auto-calculated from ``bulk_planet`` if not provided directly.
        """
        if self._stellar_dex is not None:
            return self._stellar_dex
        bulk = self.bulk_planet
        if bulk is not None:
            return conv.calculate_dex_from_bulk_planet(bulk)
        return None

    @property
    def alphas(self) -> dict[str, float] | None:
        """dict[str, float] or None : Element partitioning ratios (BSP/BP) for core formation.

        Auto-calculated from ``bulk_planet`` + ``bulk_silicate_planet`` if
        not provided directly.
        """
        if self._alphas is not None:
            return self._alphas
        if self._bulk_planet is not None and self._bulk_silicate_planet is not None:
            bp_elements = conv.convert_composition(self._bulk_planet, 'wtpt_elements')
            bsp_elements = conv.convert_composition(self._bulk_silicate_planet, 'wtpt_elements')
            return {
                el: bsp_elements[el] / bp_elements[el]
                for el in bp_elements
                if bp_elements[el] > 0 and el in bsp_elements
            }
        return None
    
    @property
    def name(self) -> str | None:
        """str or None : Planet name."""
        return self._name

    @property
    def mass(self) -> float | None:
        """float or None : Planet mass (not yet implemented)."""
        return self._mass

    @classmethod
    def from_star(cls, star: Star, **kwargs: Any) -> Planet:
        """Create a Planet from a :class:`~stellar_geology.star.Star` object.

        Parameters
        ----------
        star : :class:`~stellar_geology.star.Star`
            A Star instance with a ``stellar_dex`` composition.

        Returns
        -------
        Planet
            A new Planet initialized with the star's dex composition.
        """
        return cls(stellar_dex=star.stellar_dex, **kwargs)

    def set_alphas(self, alphas: dict[str, float] | None) -> None:
        """Update the mantle-core partitioning coefficients ``alphas``.

        Parameters
        ----------
        alphas : dict[str, float] or None
            Mantle-core partitioning coefficients, or ``None`` to clear.

        Examples
        --------
        Calculate bulk silicate planet compositions across a range of alpha values:

        >>> p = Planet.from_star(star, alphas={'Fe': 0.49, 'Ni': 0.49})
        >>> bsps = {}
        >>> for alpha_fe in (0.30, 0.40, 0.49, 0.55):
        ...     p.set_alphas({'Fe': alpha_fe, 'Ni': 0.49})
        ...     bsps[alpha_fe] = p.bulk_silicate_planet
        """
        self._alphas = alphas

    def get_composition(self, which: str, units: str = 'wtpt_oxides') -> dict[str, float]:
        """
        Return the planet's composition in the requested units.

        Output dicts always sum to the target for the requested units —
        100 for percent units (``wtpt_*``, ``molpt_*``) and 1.0 for
        fraction units (``wtfrac_*``, ``molfrac_*``).

        Parameters
        ----------
        which : str
            One of 'bulk_planet', 'bulk_silicate_planet'.
        units : str
            One of 'wtpt_oxides', 'wtpt_elements', 'wtfrac_oxides',
            'wtfrac_elements', 'molfrac_oxides', 'molfrac_elements',
            'molfrac_singleO', 'molpt_oxides', 'molpt_elements'.

        Returns
        -------
        dict[str, float]
            Composition in the requested units.

        Raises
        ------
        ValueError
            If ``which`` or ``units`` is invalid, or if the requested
            composition cannot be computed from the inputs provided to the
            Planet (e.g. ``bulk_silicate_planet`` was requested but no
            ``alphas`` were provided).

        Notes
        -----
        Planet element outputs do not include volatile elements (C, O, S).
        We assume volatiles are mostly lost during planet formation. Element
        conversions go through the oxide-based converter, which only includes
        non-volatile rock-forming elements.
        """
        valid_which = ['bulk_planet', 'bulk_silicate_planet']
        if which not in valid_which:
            raise ValueError(f"which must be one of {valid_which}, got '{which}'.")

        if units not in conv.VALID_UNITS:
            raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")

        if which == 'bulk_planet':
            base_composition = self.bulk_planet
        else:
            base_composition = self.bulk_silicate_planet

        if base_composition is None:
            raise ValueError(self._diagnose_missing_inputs(which))

        return conv.convert_composition(base_composition, units)

    def _diagnose_missing_inputs(self, which: str) -> str:
        """Build a human-readable explanation of why ``which`` cannot be
        computed from the current inputs. Used by :meth:`get_composition`
        to produce a clear ValueError message."""
        alpha_hint = ("Supply alphas via Planet(..., alphas=...), "
                      "Planet.from_star(star, alphas=...), or "
                      "planet.set_alphas(...).")
        if which == 'bulk_planet':
            if self._bulk_silicate_planet is not None and self._alphas is None:
                return ("bulk_planet cannot be computed: bulk_silicate_planet "
                        f"was provided but alphas is missing. {alpha_hint}")
            return ("bulk_planet is not set and cannot be computed. "
                    "Pass bulk_planet, stellar_dex, or "
                    "(bulk_silicate_planet + alphas) at construction.")
        # which == 'bulk_silicate_planet'
        bp_available = (self._bulk_planet is not None
                        or self._stellar_dex is not None)
        if bp_available and self._alphas is None:
            return ("bulk_silicate_planet cannot be computed: bulk_planet "
                    f"was provided but alphas is missing. {alpha_hint}")
        return ("bulk_silicate_planet is not set and cannot be computed. "
                "Pass bulk_silicate_planet, or (bulk_planet + alphas).")

    #--- CALCULATIONS BETWEEN BULK PLANET AND BULK SILICATE PLANET ---#
    def _calculate_silicate_from_bulk(self, bulk_planet: dict[str, float],
                                         alphas: dict[str, float]) -> dict[str, float]:
        """
        Calculates the bulk silicate composition given known bulk composition and
        the alphas ratio for partitioning bulk Fe between the core and mantle, as
        defined in Putirka and Rarick (2019).

        The algorithm follows the Putirka & Rarick (2019) supplementary
        spreadsheet.  Only Fe and Ni alphas are used to set their BSP
        concentrations directly.  All other elements (including Si, even if
        an alpha is supplied) are proportionally rescaled to fill the
        remaining mass (100 − Fe_BSP − Ni_BSP).  The result is converted
        from element wt% to oxide wt% and normalized to 100.

        Parameters
        ----------
        bulk_planet:    dict[str, float]
            Bulk planet composition in wt% oxides. Must contain value for FeO.
        alphas:    dict[str, float]
                Ratio of element in the bulk silicate planet and bulk planet, defined
                in Putirka and Rarick (2019): e.g., alphas = FeBSP/FeBP. Will always
                be a positive fraction <1. Used for defining which elements partition
                into a metallic core. Commonly, Fe, Si, and Ni. Fe is required when
                passing this argument: {'Fe': 0.49}.

        Returns
        -------
        dict[str, float]
            Bulk silicate planet composition in wt% oxides.
        """
        # TODO consider failure cases for other lack of keys (Ni?), units, etc...
        if "FeO" not in list(bulk_planet.keys()):
            raise ValueError("Bulk planet composition must have FeO concentration.")

        # Translate bulk_planet to wt% element basis. Reachable only via the
        # bulk_silicate_planet property after it has confirmed bulk_planet is
        # available, so this call is guaranteed to succeed.
        bulk_wtpt_elements = self.get_composition(
            which="bulk_planet", units="wtpt_elements"
        )

        # Validate alpha values
        for k, v in alphas.items():
            if v <= 0 or v >= 1:
                raise ValueError(
                    f"{k} alpha value must be a float where 0 < alpha < 1"
                )

        # --- Putirka & Rarick (2019) algorithm ---
        # 1. Fe and Ni BSP concentrations set directly from their alphas.
        fe_bsp = alphas["Fe"] * bulk_wtpt_elements["Fe"]
        ni_bsp = alphas.get("Ni", 1.0) * bulk_wtpt_elements.get("Ni", 0.0)

        # 2. Remaining mass budget for all other (lithophile) elements.
        remaining_mass = 100.0 - fe_bsp - ni_bsp

        # 3. Sum of BP element wt% for all elements except Fe and Ni.
        sum_lithophile_bp = sum(
            v for k, v in bulk_wtpt_elements.items() if k not in ("Fe", "Ni")
        )

        # 4. Build BSP element wt%: Fe/Ni are set directly; everything else
        #    is rescaled proportionally to fill the remaining mass.
        bsp_elements = {}
        for k, v in bulk_wtpt_elements.items():
            if k == "Fe":
                bsp_elements[k] = fe_bsp
            elif k == "Ni":
                bsp_elements[k] = ni_bsp
            else:
                bsp_elements[k] = remaining_mass * v / sum_lithophile_bp

        # 5. Convert BSP element wt% → wt% oxides (normalizes to 100).
        return conv.convert_to_wtpt_oxides(bsp_elements, "wtpt_elements")

    def _calculate_bulk_from_silicate(self, bulk_silicate_planet: dict[str, float],
                                      alphas: dict[str, float]) -> dict[str, float]:
        """
        Recovers a bulk planet composition from a bulk silicate planet
        composition and core partitioning fractions (alphas). This is the
        reverse of :meth:`_calculate_silicate_from_bulk`.

        The algorithm reverses the Putirka & Rarick (2019) procedure:
        Fe and Ni concentrations are inflated by their reciprocal alphas,
        and all lithophile elements are rescaled to fill the remaining
        mass budget. The result is in wt% oxides, normalized to 100.

        Note that this reversal assumes the same alpha values used in the
        forward calculation. Different alphas will produce a different
        (but internally consistent) bulk planet. This is an intrinsic
        ambiguity of the core partitioning model — see CAVEATS.md.

        Parameters
        ----------
        bulk_silicate_planet : dict[str, float]
            Bulk silicate planet composition in wt% oxides.
        alphas : dict[str, float]
            Core partitioning fractions (same definition as forward:
            alpha = BSP_element / BP_element). Must include 'Fe'.

        Returns
        -------
        dict[str, float]
            Bulk planet composition in wt% oxides.
        """
        if "Fe" not in alphas:
            raise ValueError("alphas must include 'Fe'.")

        # BSP → wt% elements
        bsp_elements = conv.convert_composition(bulk_silicate_planet, 'wtpt_elements')

        # Reverse the alpha scaling for Fe and Ni
        bp_fe = bsp_elements.get("Fe", 0.0) / alphas["Fe"]
        bp_ni = bsp_elements.get("Ni", 0.0) / alphas.get("Ni", 1.0)

        bsp_fe = bsp_elements.get("Fe", 0.0)
        bsp_ni = bsp_elements.get("Ni", 0.0)

        # Reverse the lithophile rescaling
        remaining_mass_bsp = 100.0 - bsp_fe - bsp_ni
        remaining_mass_bp = 100.0 - bp_fe - bp_ni

        if remaining_mass_bsp == 0:
            raise ValueError("BSP has no lithophile elements (Fe + Ni = 100%).")

        ratio = remaining_mass_bp / remaining_mass_bsp

        bp_elements = {}
        for el, val in bsp_elements.items():
            if el == "Fe":
                bp_elements[el] = bp_fe
            elif el == "Ni":
                bp_elements[el] = bp_ni
            else:
                bp_elements[el] = val * ratio

        return conv.convert_to_wtpt_oxides(bp_elements, "wtpt_elements")
