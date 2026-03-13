"""The Planet module defines the :class:`Planet` class for representing
planetary bulk and silicate compositions derived from stellar data or
direct geochemical inputs.
"""

from . import star
from . import conversions as conv
from . import constants as const
import warnings as w


class Planet(object):
    def __init__(self, bulk_planet=None, bulk_silicate_planet=None,
                 stellar_dex=None, alphas=None, name=None, mass=None,
                 mineralogy=None, units='wtpt_oxides'):
        """
        Returns a Planet() object.

        Parameters
        ----------
        None required. Would return an empty Planet() object with no attributes.
        Caution: do not pass multiple conflicting compositional parameters or
        it will raise an Error. Just pass one, and the others will be auto-
        matically computed for you. Like magic.

        bulk_planet : dict
            Bulk planet composition. Units specified by the `units` parameter.
        bulk_silicate_planet : dict
            Bulk silicate planet composition. Units specified by the `units`
            parameter.
        stellar_dex : dict
            Star composition in dex notation.
        alphas : dict[str: float]
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
        
        if bulk_planet is not None:
            unrecognized_keys = []
            for k in bulk_planet.keys():
                if k not in const.elements_to_oxides.keys() and k not in const.oxides_to_elements.keys():
                    unrecognized_keys.append(k)
            if len(unrecognized_keys) > 0:
                w.warn(f"{unrecognized_keys} were not recognized as compositional parameters and "
                       "will be ignored in calculations.", category=UserWarning)
        
        if bulk_planet is not None and stellar_dex is not None:
            raise ValueError("Can not pass both bulk_planet and stellar_dex.")
        
        if bulk_planet is not None and bulk_silicate_planet is not None and alphas is not None:
            raise ValueError("Cannot pass all bulk_planet, bulk_silicate_planet"
                             ", and alphas as values may be contradictory.")
        
        if bulk_silicate_planet is not None and alphas is not None and stellar_dex is not None:
            raise ValueError("Cannot pass all bulk_silicate_planet, alphas, "
                             "and stellar_dex, as values may be contradictory.")
        
    @property
    def bulk_planet(self):
        """dict or None : Bulk planet composition in wt% oxides.

        Auto-calculated from ``stellar_dex`` or
        ``bulk_silicate_planet`` + ``alphas`` if not provided directly.
        """
        if self._bulk_planet is not None:
            return self._bulk_planet
        if self._stellar_dex is not None:
            return conv.calculate_bulk_planet_from_dex(self._stellar_dex)
        if self._bulk_silicate_planet is not None and self.alphas is not None:
            return conv.calculate_bulk_from_silicate(self._bulk_silicate_planet,
                                                     self.alphas)
        # Explain what's missing
        if self._bulk_silicate_planet is not None and self._alphas is None:
            w.warn("bulk_planet cannot be computed: "
                   "bulk_silicate_planet was provided but alphas is missing.",
                   category=UserWarning)
        elif self._stellar_dex is None:
            w.warn("bulk_planet is not set and cannot be computed. "
                   "Pass bulk_planet, stellar_dex, or "
                   "(bulk_silicate_planet + alphas).",
                   category=UserWarning)
        return None

    @property
    def bulk_silicate_planet(self):
        """dict or None : Bulk silicate planet composition in wt% oxides.

        Auto-calculated from ``bulk_planet`` + ``alphas`` if not provided
        directly.
        """
        if self._bulk_silicate_planet is not None:
            return self._bulk_silicate_planet
        if self._bulk_planet is not None and self._alphas is not None:
            return self._calculate_silicate_from_bulk(bulk_planet=self._bulk_planet,
                                                      alphas=self._alphas)
        # Explain what's missing
        if self._bulk_planet is not None and self._alphas is None:
            w.warn("bulk_silicate_planet cannot be computed: "
                   "bulk_planet was provided but alphas is missing.",
                   category=UserWarning)
        elif self._bulk_planet is None:
            w.warn("bulk_silicate_planet is not set and cannot be computed. "
                   "Pass bulk_silicate_planet, or (bulk_planet + alphas).",
                   category=UserWarning)
        return None

    @property
    def stellar_dex(self):
        """dict or None : Stellar composition in dex notation.

        Auto-calculated from ``bulk_planet`` if not provided directly.
        """
        if self._stellar_dex is not None:
            return self._stellar_dex
        if self._bulk_planet is not None:
            return self._calculate_dex_from_bulk()
        if self._bulk_silicate_planet is not None and self._alphas is not None:
            self._bulk_planet = self._calculate_bulk_from_silicate()
            return self._calculate_dex_from_bulk()
        w.warn("stellar_dex is not set and cannot be computed. "
               "Pass stellar_dex, bulk_planet, or "
               "(bulk_silicate_planet + alphas).",
               category=UserWarning)
        return None

    @property
    def alphas(self):
        """dict or None : Element partitioning ratios (BSP/BP) for core formation.

        Auto-calculated from ``bulk_planet`` + ``bulk_silicate_planet`` if
        not provided directly.
        """
        if self._alphas is not None:
            return self._alphas
        if self._bulk_planet is not None and self._bulk_silicate_planet is not None:
            return self._calculate_alphas_from_bulk_and_silicate()
        if self.stellar_dex is not None and self._bulk_silicate_planet is not None:
            self._stellar_dex = self._calculate_dex_from_bulk()
            return self._calculate_alphas_from_bulk_and_silicate()
        # Explain what's missing
        missing = []
        if self._bulk_planet is None:
            missing.append("bulk_planet")
        if self._bulk_silicate_planet is None:
            missing.append("bulk_silicate_planet")
        w.warn(f"alphas cannot be computed: "
               f"{' and '.join(missing)} missing. Pass alphas directly, "
               f"or provide both bulk_planet and bulk_silicate_planet.",
               category=UserWarning)
        return None
    
    @property
    def name(self):
        """str or None : Planet name."""
        return self._name

    @property
    def mass(self):
        """float or None : Planet mass (not yet implemented)."""
        return self._mass

    @classmethod
    def from_star(cls, star):
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
        return cls(stellar_dex=star.stellar_dex)

    def get_composition(self, which, units='wtpt_oxides', normalization=None):
        """
        Return the planet's composition in the requested units with optional
        normalization.

        Parameters
        ----------
        which : str
            One of 'bulk_planet', 'bulk_silicate_planet'.
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
            Composition in the requested units, or None if the base composition
            is not available.

        Notes
        -----
        Planet element outputs do NOT include volatile elements (C, O, S).
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
        elif which == 'bulk_silicate_planet':
            base_composition = self.bulk_silicate_planet

        if base_composition is None:
            return None

        result = conv.convert_composition(base_composition, units)

        if normalization is not None and normalization != 'none':
            if units != 'molfrac_singleO':
                result = conv.normalize_composition(result, normalization, units)

        return result

    #--- CALCULATIONS BETWEEN BULK PLANET AND BULK SILICATE PLANET ---#
    def _calculate_silicate_from_bulk(self, bulk_planet, alphas):
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
        bulk_planet:    dict
            Bulk planet composition in wt% oxides. Must contain value for FeO.
        alphas:    dict{str: float}
                Ratio of element in the bulk silicate planet and bulk planet, defined
                in Putirka and Rarick (2019): e.g., alphas = FeBSP/FeBP. Will always
                be a positive fraction <1. Used for defining which elements partition
                into a metallic core. Commonly, Fe, Si, and Ni. Fe is required when
                passing this argument: {'Fe': 0.49}.

        Returns
        -------
        dict
            Bulk silicate planet composition in wt% oxides.
        """
        # TODO consider failure cases for other lack of keys (Ni?), units, etc...
        if "FeO" not in list(bulk_planet.keys()):
            raise ValueError("Bulk planet composition must have FeO concentration.")

        if self._bulk_silicate_planet is not None:
            w.warn("Warning: this Planet's bulk_silicate_planet composition "
                   "will be overwritten.", category=UserWarning)

        # Translate bulk_planet to wt% element basis
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
        self._bulk_silicate_planet = conv.convert_to_wtpt_oxides(
            bsp_elements, "wtpt_elements"
        )
        return self._bulk_silicate_planet
        
        