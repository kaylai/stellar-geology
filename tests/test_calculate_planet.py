import pytest
from stellar_geology.planet import Planet
from stellar_geology.star import Star

# Tests that the math is correct
# TODO !! add tests for intermediate calcs from stellar_dex to bulk_planet

def test_bulk_planet_from_stellar_dex():
    # known dex composition of 32768 from the Hypatia catalog (Hinkel et al.
    # 2014) as reported in Putirka and Rarick (2019)
    star_32768_dex = {
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
            'C':  -0.14,
            'O':  -0.06,
        }
    expected_bulk_planet = {
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
    p = Planet(stellar_dex=star_32768_dex)
    assert p.bulk_planet == pytest.approx(expected_bulk_planet)

def test_create_planet_from_Star():
    # known dex composition of 32768 from the Hypatia catalog (Hinkel et al.
    # 2014) as reported in Putirka and Rarick (2019)
    star_32768_dex = {
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
            'C':  -0.14,
            'O':  -0.06,
        }
    expected_bulk_planet = {
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
    s = Star(stellar_dex=star_32768_dex)
    p = Planet.from_star(s)
    
    assert p.bulk_planet == pytest.approx(expected_bulk_planet)

def test_silicate_from_bulk():
    # known bulk silicate  composition of C #5 Ringwood from the Hypatia catalog
    # (Hinkel et al. 2014) as reported in Putirka and Rarick (2019)
    c5_ringwood_wtpt_oxides = {
            'SiO2'  : 34.3228620623356,
            'TiO2'  : 0.1289398301558,
            'Al2O3' : 2.6516473296260,
            'FeO'   : 33.7355431278904,
            'MgO'   : 23.9923487223086,
            'CaO'   : 2.0517814044420,
            'Na2O'  : 0.7781373566196,
            'Cr2O3' : 0.4498165616991,
            'NiO'   : 1.8889236049229,									
            }
    # alphas from Putirka and Rarick (2019) Supplementary spreadsheet 2
    alphas = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}
    expected_silicate_planet = {
        'SiO2'  : 45.0623584204862,
        'TiO2'  : 0.1692846252333,
        'Al2O3' : 3.4813379535577,
        'FeO'   : 15.3424822627129,
        'MgO'   : 31.4994657354165,
        'CaO'   : 2.6937762031481,
        'Na2O'  : 1.0216136521680,
        'Cr2O3' : 0.5905624970884,
        'NiO'   : 0.1391186501890,
    }
    p = Planet(bulk_planet=c5_ringwood_wtpt_oxides, alphas=alphas)
    # rel=1e-4 accounts for small differences in molecular weight constants
    # between our code and the Putirka & Rarick (2019) spreadsheet
    assert p.bulk_silicate_planet == pytest.approx(expected_silicate_planet, rel=1e-4)

def test_mantle_CIPW_norm():
    # TODO test for calculating CIPW norm from mantle (bulk_silicate_planet) composition.
    pass


# ============================================================================
# REVERSE PIPELINE: BSP + alphas → bulk_planet
# ============================================================================
def test_bulk_planet_from_bsp_and_alphas():
    """BSP + alphas should recover bulk_planet via the reverse calculation."""
    c5_ringwood_wtpt_oxides = {
        'SiO2'  : 34.3228620623356,
        'TiO2'  : 0.1289398301558,
        'Al2O3' : 2.6516473296260,
        'FeO'   : 33.7355431278904,
        'MgO'   : 23.9923487223086,
        'CaO'   : 2.0517814044420,
        'Na2O'  : 0.7781373566196,
        'Cr2O3' : 0.4498165616991,
        'NiO'   : 1.8889236049229,
    }
    alphas = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}

    # Forward: bulk → BSP
    p_fwd = Planet(bulk_planet=c5_ringwood_wtpt_oxides, alphas=alphas)
    bsp = p_fwd.bulk_silicate_planet

    # Reverse: BSP + alphas → bulk
    p_rev = Planet(bulk_silicate_planet=bsp, alphas=alphas)
    recovered_bulk = p_rev.bulk_planet

    assert recovered_bulk == pytest.approx(c5_ringwood_wtpt_oxides, rel=1e-4)


def test_stellar_dex_from_bulk_planet():
    """Planet initialized with bulk_planet should compute stellar_dex in reverse."""
    star_32768_dex = {
        'Si': 0.27,
        'Ti': 4.61436e-16,
        'Cr': 0.08,
        'Al': 0.23,
        'Fe': 0.02,
        'Mn': 0.06,
        'Mg': 0.21,
        'Ca': 0.1,
        'Na': 0.3,
        'Ni': 0.04,
        'C':  -0.14,
        'O':  -0.06,
    }
    p = Planet(stellar_dex=star_32768_dex)
    bulk = p.bulk_planet

    # New planet from bulk, should compute dex in reverse
    p2 = Planet(bulk_planet=bulk)
    recovered_dex = p2.stellar_dex

    # Relative dex differences should be preserved
    common = [el for el in star_32768_dex
              if el in recovered_dex and el not in ('C', 'O', 'S')]
    ref = 'Fe'
    for el in common:
        if el == ref:
            continue
        original_diff = star_32768_dex[el] - star_32768_dex[ref]
        recovered_diff = recovered_dex[el] - recovered_dex[ref]
        assert recovered_diff == pytest.approx(original_diff, abs=1e-6)


def test_alphas_from_bulk_and_bsp():
    """Planet with both bulk and BSP should compute alphas."""
    c5_ringwood = {
        'SiO2'  : 34.3228620623356,
        'TiO2'  : 0.1289398301558,
        'Al2O3' : 2.6516473296260,
        'FeO'   : 33.7355431278904,
        'MgO'   : 23.9923487223086,
        'CaO'   : 2.0517814044420,
        'Na2O'  : 0.7781373566196,
        'Cr2O3' : 0.4498165616991,
        'NiO'   : 1.8889236049229,
    }
    alphas_orig = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}

    p_fwd = Planet(bulk_planet=c5_ringwood, alphas=alphas_orig)
    bsp = p_fwd.bulk_silicate_planet

    # Create planet with both bulk and BSP, no alphas
    p_rev = Planet(bulk_planet=c5_ringwood, bulk_silicate_planet=bsp)
    computed_alphas = p_rev.alphas

    assert computed_alphas["Fe"] == pytest.approx(alphas_orig["Fe"], rel=1e-4)
    assert computed_alphas["Ni"] == pytest.approx(alphas_orig["Ni"], rel=1e-4)