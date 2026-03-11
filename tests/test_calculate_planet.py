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
    C5Ringwood_wtptOxides = {
            'SiO2'  : 34.32286206,
            'TiO2'  : 0.12893983,
            'Cr2O3' : 2.65164733,
            'Al2O3' : 33.73554313,
            'FeO'   : 23.99234872,
            'MgO'   : 2.051781404,
            'NiO'   : 0.778137357,
            'CaO'   : 0.449816562,
            'Na2O'  : 1.888923605,									
            }
    # alphas from Putirka and Rarick (2019) Supplementary spreadsheet 2
    alphas = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}
    expected_silicate_planet = {
        'SiO2'  : 45.06235842,
        'TiO2'  :0.169284625,
        'Al2O3' :3.481337954,
        'FeO'  : 15.34248226,
        'MgO'  : 31.49946574,
        'CaO'  : 2.693776203,
        'Na2O'  : 1.021613652,
        'Cr2O3'  : 0.590562497,
        'NiO'  : 0.13911865,
    }
    p = Planet(bulk_planet=C5Ringwood_wtptOxides, alphas=alphas)
    assert p.bulk_silicate_planet == pytest.approx(expected_silicate_planet)
