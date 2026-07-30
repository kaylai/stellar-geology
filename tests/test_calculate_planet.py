import pytest
from stellar_geology.planet import Planet
from stellar_geology.star import Star
from stellar_geology import constants as cnst

# Tests that the math is correct
# TODO !! add tests for intermediate calcs from stellar_dex to bulk_planet

# known dex composition of 32768 from the Hypatia catalog (Hinkel et al.
# 2014) as reported in Putirka and Rarick (2019)
STELLAR_DEX_32768 = {
    'Si': 0.27,
    'Ti': 4.61436*10**(-16),
    'Cr': 0.08,
    'Al': 0.23,
    'Fe': 0.02,
    'Mn': 0.06,
    'Mg': 0.21,
    'Ca': 0.1,
    'Na': 0.3,
    'Ni': 0.04,
    }

BULK_PLANET = {
    'SiO2' : 42.37190416490457,
    'TiO2' : 0.07255698934653063,
    'Cr2O3': 0.44576451274180956,
    'Al2O3': 2.727237231195202,
    'FeO'  : 24.250084830232957,
    'MnO'  : 0.28131691827503635,
    'MgO'  : 25.332083338986063,
    'NiO'  : 1.484564845659013,
    'CaO'  : 1.68697066143704,
    'Na2O' : 1.347516507221789
    }

# alphas from Putirka and Rarick (2019) Supplementary spreadsheet 2 for c5_ringwood
ALPHAS = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}

BULK_SILICATE_PLANET = {
    'SiO2': 40.10963920477409,
    'FeO': 11.57137343438085,
    'NiO': 0.11471860511044625,
    'TiO2': 0.10966380356353714,
    'Cr2O3': 0.673735671796463,
    'Al2O3': 4.121990323558761,
    'MnO': 0.4251869260654258,
    'MgO': 38.28731919779622,
    'CaO': 2.5497146573947225,
    'Na2O': 2.036658175559482
    }

def test_bulk_planet_from_stellar_dex():
    expected_bulk_planet = BULK_PLANET
    p = Planet(stellar_dex=STELLAR_DEX_32768)
    assert p.bulk_planet == pytest.approx(expected_bulk_planet)

def test_create_planet_from_Star():
    # known dex composition of 32768 from the Hypatia catalog (Hinkel et al.
    # 2014) as reported in Putirka and Rarick (2019)
    expected_bulk_planet = BULK_PLANET
    s = Star(stellar_dex=STELLAR_DEX_32768)
    p = Planet.from_star(s)
    assert p.bulk_planet == pytest.approx(expected_bulk_planet)

def test_silicate_from_bulk():    
    p = Planet(bulk_planet=BULK_PLANET, alphas=ALPHAS)
    # rel=1e-4 accounts for small differences in molecular weight constants
    # between our code and the Putirka & Rarick (2019) spreadsheet
    assert p.bulk_silicate_planet == pytest.approx(BULK_SILICATE_PLANET, rel=1e-4)

def test_mantle_CIPW_norm():
    # TODO test for calculating CIPW norm from mantle (bulk_silicate_planet) composition.
    pass


# ============================================================================
# REVERSE PIPELINE: BSP + alphas → bulk_planet
# ============================================================================
def test_bulk_planet_from_bsp_and_alphas():
    """BSP + alphas should recover bulk_planet via the reverse calculation."""  
    # Forward: bulk to BSP
    p_fwd = Planet(bulk_planet=BULK_PLANET, alphas=ALPHAS)
    bsp = p_fwd.bulk_silicate_planet

    # Reverse: BSP + alphas → bulk
    p_rev = Planet(bulk_silicate_planet=bsp, alphas=ALPHAS)
    recovered_bulk = p_rev.bulk_planet

    assert recovered_bulk == pytest.approx(BULK_PLANET, rel=1e-4)

def test_stellar_dex_from_bulk_planet():
    """Planet initialized with bulk_planet should compute stellar_dex in reverse."""
    p = Planet(stellar_dex=STELLAR_DEX_32768)
    bulk = p.bulk_planet

    # New planet from bulk, should compute dex in reverse
    p2 = Planet(bulk_planet=bulk)
    recovered_dex = p2.stellar_dex

    # Relative dex differences should be preserved
    common = [el for el in STELLAR_DEX_32768
              if el in recovered_dex and el not in ('C', 'O', 'S')]
    ref = 'Fe'
    for el in common:
        if el == ref:
            continue
        original_diff = STELLAR_DEX_32768[el] - STELLAR_DEX_32768[ref]
        recovered_diff = recovered_dex[el] - recovered_dex[ref]
        assert recovered_diff == pytest.approx(original_diff, abs=1e-6)

def test_alphas_from_bulk_and_bsp():
    """Planet with both bulk and BSP should compute alphas."""
    p_fwd = Planet(bulk_planet=BULK_PLANET, alphas=ALPHAS)
    bsp = p_fwd.get_composition(which="bulk_silicate_planet")

    # Create planet with both bulk and BSP, no alphas
    p_rev = Planet(bulk_planet=BULK_PLANET, bulk_silicate_planet=bsp)
    computed_alphas = p_rev.alphas

    for k, v in ALPHAS.items():
        assert v == pytest.approx(computed_alphas[k], rel=1e-6)

def test_all_alphas_used():
    """Ensure that any values in the alphas dict are not ignored."""
    c5_ringwood = {
        'SiO2'  : 34.315267345094774,
        'TiO2'  : 0.15529780644418686,
        'Al2O3' : 2.6496934995271215,
        'FeO'   : 33.73568069951041,
        'MgO'   : 23.985488843225028,
        'CaO'   : 2.0490120369450553,
        'Na2O'  : 0.7778814840992446,
        'Cr2O3' : 0.4443874843135551,
        'NiO'   : 1.887290800840637,
    }
    
    expected_bsp_wtpt_elements = {
        'Si'    : 10.5829853,
        'Ti'    : 0.09211496,
        'Al'    : 1.85048719,
        'Fe'    : 43.253091,
        'Mg'    : 28.6293294,
        'Ca'    : 3.3816251,
        'Na'    : 1.52296733,
        'Cr'    : 0.9027266,
        'Ni'    : 9.78465569,   
    }
    
    # create unique fake alphas for all known cations using their mass
    set_alphas = {
        'Si'    : 0.1,
        'Ti'    : 0.15,
        'Al'    : 0.2,
        'Fe'    : 0.25,
        'Mg'    : 0.3,
        'Ca'    : 0.35,
        'Na'    : 0.4,
        'Cr'    : 0.45,
        'Ni'    : 1,
    }        
    
    p = Planet(bulk_planet=c5_ringwood, alphas=set_alphas)
    bsp = p.get_composition(which="bulk_silicate_planet", units="wtpt_elements")
    
    assert bsp == pytest.approx(expected_bsp_wtpt_elements, rel=1e-4)

def test_bad_alpha_raises():
    """Ensure that any key in the alphas dict that is not a known cation raises a
    ValueError. Ensure that any value in the alphas dict is a float >0 and <=1"""
    # known dex composition of 32768 from the Hypatia catalog (Hinkel et al.
    # 2014) as reported in Putirka and Rarick (2019)
    star = Star(STELLAR_DEX_32768)
    alphas = {"Fe": 0.4, "Ni": 0.2, "fakeelement": 0.1}
    alphas_bad_Ni = {"Fe": 0.4, "Ni": 0}
    
    with pytest.raises(ValueError, match="Alpha values contain"):
        Planet.from_star(star, alphas=alphas)
    with pytest.raises(ValueError, match="Ni alpha value must"):
        Planet.from_star(star, alphas=alphas_bad_Ni)

def test_zero_value_in_comp_returns_None_alpha():
    """If alphas are calculated from a BSP plus either a stellar_dex or BP, any compositional
    component with a 0 value in the BSP, stellar_dex, or BP will result in an alpha value of 0
    for that element. check_alphas() would then raise since alpha cannot be =0. Instead, the
    alphas @property in planet.py pop's any key <=0 from the derived alphas dict. This behavior
    should only apply to lazy computing alpha values from Planet.alphas. If a user passes alpha
    values with a 0 value in the dict, it should raise as in `test_bad_alpha_raises()` in this test
    file."""
    bsp_w_zero = BULK_SILICATE_PLANET
    bsp_w_zero['P2O5'] = 0
    
    astar = Star(STELLAR_DEX_32768)
    p = Planet.from_star(astar, bulk_silicate_planet=bsp_w_zero)
    computed_alphas = p.alphas
    print(computed_alphas)
    
    assert 'P' not in computed_alphas.keys()
    for k in computed_alphas.keys():
        if k not in ALPHAS:
            assert computed_alphas[k] == 1
        else:
            assert computed_alphas[k] == pytest.approx(ALPHAS[k], rel=1e-6)
    
