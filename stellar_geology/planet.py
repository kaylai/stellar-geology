from . import star
from . import conversions as conv
from . import constants as const

class Planet(object):
    def __init__(self, bulk_planet=None, bulk_silicate_planet=None, 
                 stellar_dex=None, alphaFe=None, name=None, mass=None):
        """
        Returns a Planet() object.
        
        Parameters
        ----------
        None required. Would return an empty Planet() object with no attributes.
        Caution: do not pass multiple conflicting compositional parameters or
        it will raise an Error. Just pass one, and the others will be auto-
        matically computed for you. Like magic.
        
        bulk_planet:    dict
            Bulk planet composition in wt% oxides.
        bulk_silicate_planet: dict
            Bulk silicate planet composition in wt% oxides.
        stellar_dex:    dict
            Star composition in dex notation.
        alphaFe:    float
            Ratio of Fe in the bulk silicate planet and bulk planet, defined
            in Putirka and Rarick (2019): alphaFe = FeBSP/FeBP. Will always
            be a positive fraction <1.
        name: str
            Arbitrary name for your planet as a string. Can be anything. Either
            aa, either bb, even zombocom. At zombocom you can do anything.
        mass:   float
            Planet mass in some units that I don't know because this isn't
            implemented yet. So, put whatever float you want here. It won't
            make any difference. 
        
        Returns
        -------
        Planet() object.
        """
        self._bulk_planet = bulk_planet
        self._bulk_silicate_planet = bulk_silicate_planet
        self._stellar_dex = stellar_dex
        self._alphaFe = alphaFe
        self._name = name
        self._mass = mass
        
        if bulk_planet is not None:
            for k in bulk_planet.keys():
                if k not in list(const.oxides_to_elements.keys()):
                    raise ValueError("bulk_planet must be passed as oxides (e.g., 'SiO2 not 'Si').")
        
        if bulk_planet is not None and stellar_dex is not None:
            raise ValueError("Can not pass both bulk_planet and stellar_dex.")
        
        if bulk_planet is not None and bulk_silicate_planet is not None and alphaFe is not None:
            raise ValueError("Cannot pass all bulk_planet, bulk_silicate_planet"
                             ", and alphaFe as values may be contradictory.")
        
        if bulk_silicate_planet is not None and alphaFe is not None and stellar_dex is not None:
            raise ValueError("Cannot pass all bulk_silicate_planet, alphaFe, "
                             "and stellar_dex, as values may be contradictory.")
        
    @property
    def bulk_planet(self):
        if self._bulk_planet is not None:
            return self._bulk_planet
        if self._stellar_dex is not None:
            return conv.calculate_bulk_planet_from_dex(self._stellar_dex)
        if self._bulk_silicate_planet is not None and self.alphaFe is not None:
            return conv.calculate_bulk_from_silicate(self._bulk_silicate_planet,
                                                     self.alphaFe)
        return None
    
    @property
    def bulk_silicate_planet(self):
        if self._bulk_silicate_planet is not None:
            return self._bulk_silicate_planet
        if self._bulk_planet is not None and self._alphaFe is not None:
            return self._calculate_silicate_from_bulk()
        return None
    
    @property
    def stellar_dex(self):
        if self._stellar_dex is not None:
            return self._stellar_dex
        if self._bulk_planet is not None:
            return self._calculate_dex_from_bulk()
        if self._bulk_silicate_planet is not None and self._alphaFe is not None:
            self._bulk_planet = self._calculate_bulk_from_silicate()
            return self._calculate_dex_from_bulk()
        return None
    
    @property
    def alphaFe(self):
        if self._alphaFe is not None:
            return self._alphaFe
        if self._bulk_planet is not None and self._bulk_silicate_planet is not None:
            return self._calculate_alphaFe_from_bulk_and_silicate()
        if self.stellar_dex is not None and self._bulk_silicate_planet is not None:
            self._stellar_dex = self._calculate_dex_from_bulk()
            return self._calculate_alphaFe_from_bulk_and_silicate()
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
        
        
        
        