import pytest
from stellar_geology.planet import Planet

# Tests that function works as expected logistically
def test_empty_planet_returns_none():
    p = Planet()
    assert p.bulk_planet is None
    assert p.bulk_silicate_planet is None
    assert p.stellar_dex is None
    assert p.alphaFe is None
    assert p.name is None
    assert p.mass is None

def test_bulk_planet_returns_directly():
    comp = {"SiO2": 45.0, "MgO": 30.0}
    p = Planet(bulk_planet=comp)
    assert p.bulk_planet == comp

def test_name_returns_directly():
    p = Planet(name="zombocom")
    assert p.name == "zombocom"

def test_conflicting_inputs_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45}, stellar_dex={"Si": 0.1})

def test_all_three_bulk_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45.0},
               bulk_silicate_planet={"SiO2": 50.0},
               alphaFe=0.5)

def test_bsp_alpha_dex_raises():
    with pytest.raises(ValueError):
        Planet(bulk_silicate_planet={"SiO2": 50.0},
               alphaFe=0.5,
               stellar_dex={"Si": 0.1})

# Tests that the math is correct
# TODO this will be moved to a test_conversions.py file eventually...
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
    p = Planet(stellar_dex=star_32768_dex)
    expected_bulk_planet = {'SiO2': 42.37190416490457,
            'TiO2': 0.07255698934653063,
            'Cr2O3': 0.44576451274180956,
            'Al2O3': 2.727237231195202,
            'FeO': 24.250084830232957,
            'MnO': 0.28131691827503635,
            'MgO': 25.332083338986063,
            'NiO': 1.484564845659013,
            'CaO': 1.68697066143704,
            'Na2O': 1.347516507221789
            }
    assert p.bulk_planet == pytest.approx(expected_bulk_planet)
