from . import star
from . import conversions as conv
from . import constants as const

class Planet(object):
    def __init__(self, bulk_planet=None, bulk_silicate_planet=None,
                 stellar_dex=None, alphas=None, name=None, mass=None,
                 units='wtpt_oxides'):
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
            for k in bulk_planet.keys():
                if k not in list(const.oxides_to_elements.keys()):
                    raise ValueError("bulk_planet must be passed as oxides (e.g., 'SiO2 not 'Si').")
        
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
        if self._bulk_planet is not None:
            return self._bulk_planet
        if self._stellar_dex is not None:
            return conv.calculate_bulk_planet_from_dex(self._stellar_dex)
        if self._bulk_silicate_planet is not None and self.alphas is not None:
            return conv.calculate_bulk_from_silicate(self._bulk_silicate_planet,
                                                     self.alphas)
        return None
    
    @property
    def bulk_silicate_planet(self):
        if self._bulk_silicate_planet is not None:
            return self._bulk_silicate_planet
        if self._bulk_planet is not None and self._alphas is not None:
            return self._calculate_silicate_from_bulk(bulk_planet=self._bulk_planet,
                                                      alphas=self._alphas)
        return None
    
    @property
    def stellar_dex(self):
        if self._stellar_dex is not None:
            return self._stellar_dex
        if self._bulk_planet is not None:
            return self._calculate_dex_from_bulk()
        if self._bulk_silicate_planet is not None and self._alphas is not None:
            self._bulk_planet = self._calculate_bulk_from_silicate()
            return self._calculate_dex_from_bulk()
        return None
    
    @property
    def alphas(self):
        if self._alphas is not None:
            return self._alphas
        if self._bulk_planet is not None and self._bulk_silicate_planet is not None:
            return self._calculate_alphas_from_bulk_and_silicate()
        if self.stellar_dex is not None and self._bulk_silicate_planet is not None:
            self._stellar_dex = self._calculate_dex_from_bulk()
            return self._calculate_alphas_from_bulk_and_silicate()
        return None
    
    @property
    def name(self):
        return self._name
    
    @property
    def mass(self):
        return self._mass
    
    @classmethod
    def from_star(cls, star):
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
        
        # must translate bulk_planet to be on wt% cation basis
        for k, v in alphas.items():
            if v <= 0 or v >= 1:
                raise ValueError(f"{k} alpha value must be a float where 0 < alpha < 1")
        
        partitioned_silicate_concentrations = {k:v*bulk_planet[k] for k, v in alphas.items()}
        mantle_only_concentrations = {k:v for k, v in bulk_planet.items() if k not in alphas.keys()}
        silicate = partitioned_silicate_concentrations
        pass

    # def _calculate_bulk_from_silicate(bulk_silicate_planet, alphas):
    #     """
    #     Calculates the bulk planet composition given known bulk silicate composition
    #     and the alphas ratio for partitioning bulk Fe between the core and mantle,
    #     as defined in Putirka and Rarick (2019).
        
    #     Parameters
    #     ----------
    #     bulk_silicate_planet:    dict
    #         Bulk silicate planet composition in wt% oxides.
    #     alphas:    float
    #         Ratio of Fe in the bulk silicate planet and bulk planet, defined
    #         in Putirka and Rarick (2019): alphas = FeBSP/FeBP. Will always
    #         be a positive fraction <1.
        
    #     Returns
    #     -------
    #     dict
    #         Bulk planet composition in wt% oxides.
    #     """
    #     planet_C5Ringwood_silicate = planet_C5Ringwood.get_composition(which="bulk_planet", units="wtpt_elements")
    #     for k, v in alphas.items():
    #         planet_C5Ringwood_silicate[k] = v*planet_C5Ringwood_silicate[k]
    #     return alphas * bulk_planet  
        
        