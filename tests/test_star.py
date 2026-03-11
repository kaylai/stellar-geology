import pytest
from stellar_geology.star import Star

# Tests that function works as expected logistically
def test_empty_star_returns_none():
    s = Star()
    assert s.stellar_dex is None
    assert s.ax is None
    assert s.atomsRefSolar is None
    assert s.totalWtAtoms is None
    assert s.wtptElements is None
    assert s.wtptOxides is None
    assert s.name is None
    assert s.mass is None

def test_stellar_dex_returns_directly():
    comp = {"Si": 25.0, "Mg": 23.0}
    s = Star(stellar_dex=comp)
    assert s.stellar_dex == comp

def test_wrong_dex_units_raises():
    with pytest.raises(ValueError):
        Star(stellar_dex={"SiO2": 45.0})

def test_unrecognized_units_warns():
    with pytest.warns(UserWarning, match="not recognized"):
        Star(stellar_dex={"Si": 25.0, "Zq": 0.5})

def test_name_returns_directly():
    s = Star(name="zombocom")
    assert s.name == "zombocom"

def test_mass_returns_directly():
    s = Star(mass=1000.1)
    assert s.mass == 1000.1
