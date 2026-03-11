import pytest
from stellar_geology.planet import Planet

# Tests that function works as expected logistically
def test_empty_planet_returns_none():
    p = Planet()
    assert p.bulk_planet is None
    assert p.bulk_silicate_planet is None
    assert p.stellar_dex is None
    assert p.alphas is None
    assert p.name is None
    assert p.mass is None

def test_bulk_planet_returns_directly():
    comp = {"SiO2": 45.0, "MgO": 30.0}
    p = Planet(bulk_planet=comp)
    assert p.bulk_planet == comp

def test_name_returns_directly():
    p = Planet(name="zombocom")
    assert p.name == "zombocom"

def test_mass_returns_directly():
    p = Planet(mass=1000.2)
    assert p.mass == 1000.2

def test_conflicting_inputs_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45}, stellar_dex={"Si": 0.1})

def test_all_three_bulk_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45.0},
               bulk_silicate_planet={"SiO2": 50.0},
               alphas={'Fe':0.5})

def test_bsp_alpha_dex_raises():
    with pytest.raises(ValueError):
        Planet(bulk_silicate_planet={"SiO2": 50.0},
               alphas={'Fe':0.5},
               stellar_dex={"Si": 0.1})

def test_calculate_silicate_from_bulk_noFe():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45}, alphas=[]).bulk_silicate_planet