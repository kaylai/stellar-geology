"""The Planet module defines the :class:`Planet` class for representing
planetary bulk and silicate compositions derived from stellar data or
direct geochemical inputs.
"""

from __future__ import annotations

import warnings
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
        
        if alphas is not None:
            const.check_alphas(alphas)
        
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

    def _fractionate_core(self) -> tuple[dict[str, float], float, float,
                                       dict[str, float], dict[str, float]] | None:
        """Solve the core partitioning mass balance from BP and BSP.

        Computes the element-wise concentration ratios BSP/BP on the normalized wt%
        element basis and splits them into alphas (0 < r <= 1, elements partitioned into
        the core) and lithophile elements (r > 1, assumed fully retained in the
        silicate).

        Splitting the mass of the planet between the core and mantle depletes the mantle
        in some elements and enriches them in all others. Any elements not partitioned
        into the core have a common enrichment factor of 1/f, where f is the
        silicate mass fraction. That enrichment factor is used to compute core mass
        properties.

        Returns
        -------
        tuple or None
            ``(alphas, enrichment, spread, ratios, bp_elements)`` where
            ``alphas`` is the derived alpha dict, ``enrichment`` is the
            mass-weighted lithophile enrichment factor (1/f), ``spread``
            is the absolute spread of individual lithophile ratios around
            it, ``ratios`` holds r for every element with BP > 0 (r = 0
            for elements absent from the BSP), and ``bp_elements`` is the
            BP in normalized wt% elements. Returns None if either
            composition is unavailable.

        Warns
        -----
        UserWarning
            If the lithophile ratios disagree beyond our set tolerance of 1e-3
            (relative), meaning the BP/BSP pair is not exactly representable by the core
            partitioning model and all derived quantities are best-fit approximations.
        """
        bulk = self.bulk_planet
        silicate = self.bulk_silicate_planet
        if bulk is None or silicate is None:
            return None
        bp_elements = conv.convert_composition(bulk, 'wtpt_elements')
        bsp_elements = conv.convert_composition(silicate, 'wtpt_elements')
        ratios = {
            el: bsp_elements.get(el, 0.0) / bp_elements[el]
            for el in bp_elements
            if bp_elements[el] > 0
        }

        # Ratios in (0, 1] are alphas: those elements partition into the core.
        # Ratios of 0 mean the element is missing from one composition; they
        # cannot be represented as alphas (check_alphas requires alpha > 0).
        alphas = {el: r for el, r in ratios.items() if 0 < r <= 1}

        # Ratios > 1 are lithophile elements enriched by core removal. For a
        # model-consistent BP/BSP pair they all share the same enrichment
        # factor; a spread beyond tolerance means no alpha set can exactly
        # reproduce this pair.
        enrichment_spread_tolerance = 1e-3
        lithophiles = {el: r for el, r in ratios.items() if r > 1}
        if lithophiles:
            enrichment = (sum(bsp_elements[el] for el in lithophiles)
                          / sum(bp_elements[el] for el in lithophiles))
            spread = max(lithophiles.values()) - min(lithophiles.values())
            if spread / enrichment > enrichment_spread_tolerance:
                warnings.warn(
                    "Lithophile enrichment factors are inconsistent across "
                    f"elements (spread {spread:.4f} around {enrichment:.4f}). "
                    "This bulk planet / bulk silicate planet pair cannot be "
                    "exactly reproduced by any alpha set; derived alphas, "
                    "core_mass_fraction, and core_composition are best-fit "
                    "approximations.")
        else:
            # No enriched elements: BP and BSP are identical, so no core.
            enrichment = 1.0
            spread = 0.0
        return alphas, enrichment, spread, ratios, bp_elements

    @property
    def alphas(self) -> dict[str, float] | None:
        """dict[str, float] or None : Element partitioning ratios (BSP/BP) for core formation.

        Auto-calculated from ``bulk_planet`` + ``bulk_silicate_planet`` if
        not provided directly.

        Alphas are concentration ratios on the normalized wt% element
        basis, not mass fractions: the fraction of an element's mass
        residing in the silicate is alpha * silicate_mass_fraction.
        Elements absent from the dict are treated as fully lithophile.
        """
        if self._alphas is not None:
            const.check_alphas(self._alphas)
            return self._alphas
        core_frac_params = self._fractionate_core()
        if core_frac_params is None:
            return None
        calculated_alphas = core_frac_params[0]
        const.check_alphas(calculated_alphas)
        return calculated_alphas

    @property
    def core_mass_fraction(self) -> float | None:
        """float or None : Core mass fraction (CMF) of the planet.

        Derived from the lithophile enrichment factor: elements with no
        core affinity are concentrated in the silicate by exactly the
        mass removed into the core, so their common BSP/BP ratio fixes
        the core size. Requires both ``bulk_planet`` and
        ``bulk_silicate_planet`` to be available (directly or derived);
        returns None otherwise.

        Caveats: the value is a mass fraction on the volatile-free,
        oxygen-free element basis, not the true planet mass basis —
        converting to a true CMF requires accounting for oxygen bound in
        the mantle oxides and light elements in the core. It also
        inherits the anchor assumption that the most silicate-enriched
        elements are perfectly lithophile; a uniform depletion of all
        elements into the core is invisible to normalized compositions.
        See CAVEATS.md.
        """
        core_frac_params = self._fractionate_core()
        if core_frac_params is None:
            return None
        enrichment = core_frac_params[1]
        return (enrichment - 1.0) / enrichment

    @property
    def silicate_mass_fraction(self) -> float | None:
        """float or None : Silicate (mantle) mass fraction of the planet.

        The complement of ``core_mass_fraction``; see that property for
        derivation and caveats. Multiply an alpha by this value to get
        the fraction of that element's mass residing in the silicate.
        """
        cmf = self.core_mass_fraction
        if cmf is None:
            return None
        return 1.0 - cmf

    @property
    def core_composition(self) -> dict[str, float] | None:
        """dict[str, float] or None : Core composition in wt% elements (metal basis).

        Per-element core masses follow from the mass balance: the
        fraction of an element's planetary inventory in the core is
        1 - alpha / enrichment, where enrichment is the common
        lithophile BSP/BP ratio. Returns an empty dict if the planet has
        no core (BP equals BSP), or None if the compositions needed are
        unavailable. Shares the element-basis and lithophile-anchor
        caveats of ``core_mass_fraction``.
        """
        core_frac_params = self._fractionate_core()
        if core_frac_params is None:
            return None
        _, enrichment, _, ratios, bp_elements = core_frac_params
        core = {}
        for el, r in ratios.items():
            # retention within float noise of 1 means fully lithophile
            retention = min(r / enrichment, 1.0)
            if retention >= 1.0 - 1e-12:
                continue
            core[el] = (1.0 - retention) * bp_elements[el]
        total = sum(core.values())
        if total == 0:
            return {}
        return {el: 100.0 * m / total for el, m in core.items()}
    
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
        if alphas is None:
            # allows ._alpha attr to be cleared
            self._alphas = alphas
        else:
            const.check_alphas(alphas)
            self._alphas = alphas
    
    def set_bulk_silicate_planet(self, bulk_silicate_planet: dict[str, float] | None,
                                 units: str = 'wtpt_oxides') -> None:
        """Update the bulk_silicate_planet composition attribute
        
        Parameters
        ----------
        bulk_silicate_planet: dict[str, float] or None
            Bulk silicate planet composition. Units specified by the `units` parameter.
        units : str
                    Units of the bulk_planet and bulk_silicate_planet dicts. Defaults
                    to 'wtpt_oxides'. Any valid unit string is accepted (e.g.,
                    'wtpt_elements', 'molfrac_oxides'). Compositions are converted
                    to wt% oxides internally on init.
        """
        if units not in conv.VALID_UNITS:
            raise ValueError(f"units must be one of {conv.VALID_UNITS}, got '{units}'.")
        self._bulk_silicate_planet = bulk_silicate_planet

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
                in Putirka and Rarick (2019): e.g., alphas = FeBSP/FeBP on a cation wt%
                basis. Alpha will always be a positive fraction <1. Used for defining
                which elements partition into a metallic core. Commonly, Fe, Si, and Ni.
                Fe is required when passing this argument: {'Fe': 0.49}.

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
        const.check_alphas(alphas)

        # --- Putirka & Rarick (2019) algorithm ---
        # 1. BSP concentrations set directly from their alphas.        
        bsp_elements = {}
        for cation in const.cationMass.keys():
            if cation in alphas.keys():
                bsp_conc = alphas.get(cation, 1.0) * bulk_wtpt_elements.get(cation, 0.0)
                bsp_elements[cation] = bsp_conc

        # 2. Remaining mass budget for all other (lithophile) elements.
        remaining_mass = 100.0 - sum(bsp_elements.values())

        # 3. Sum of BP element wt% for all elements not partitioned into core.
        sum_lithophile_bp = sum(
            v for k, v in bulk_wtpt_elements.items() if k not in bsp_elements.keys()
        )

        # 4. Build BSP element wt%: Any elements w/alphas are set directly; everything
        # else is rescaled proportionally to fill the remaining mass.
        for k, v in bulk_wtpt_elements.items():
            if k not in bsp_elements.keys():
                bsp_elements[k] = remaining_mass * v / sum_lithophile_bp

        # 5. Convert BSP wt% elements to wt% oxides (normalizes to 100).
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
        if alphas is None:
            alphas = {k: 1 for k in const.cationMass.keys()}
            
        if "Fe" not in alphas:
            raise ValueError("alphas must include 'Fe'.")

        # BSP wt% oxides to wt% elements
        bsp_elements = conv.convert_composition(bulk_silicate_planet, 'wtpt_elements')

        # Reverse the alpha scaling for any elements partitioned into core
        bp_elements = {}
        for cation in const.cationMass.keys():
            if cation in alphas.keys():
                bp_conc = bsp_elements.get(cation, 0.0) / alphas.get(cation, 1.0)
                bp_elements[cation] = bp_conc

        # Reverse the lithophile rescaling
        remaining_mass  = 100.0 - sum(bp_elements.values())
        
        sum_siderophile_bsp = sum(
            v for k, v in bsp_elements.items() if k not in bp_elements.keys()
        )

        # if no alphas are passed or none <1, just return bulk composition
        if remaining_mass == 0:
            return conv.convert_to_wtpt_oxides(bp_elements, "wtpt_elements")
        
        # 4. Build BP element wt%: Any elements w/alphas are set directly; everything
        # else is rescaled proportionally to fill the remaining mass.
        for k, v in bsp_elements.items():
            if k not in bp_elements.keys():
                bp_elements[k] = remaining_mass * v / sum_siderophile_bsp

        return conv.convert_to_wtpt_oxides(bp_elements, "wtpt_elements")
