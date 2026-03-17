import pytest
from stellar_geology.star import Star

# Tests that function works as expected logistically
def test_empty_star_returns_none():
    s = Star()
    assert s.stellar_dex is None
    assert s.ax is None
    assert s.atoms_ref_solar is None
    assert s.total_wt_atoms is None
    assert s.wtpt_elements is None
    assert s.wtpt_oxides is None
    assert s.name is None
    assert s.mass is None

def test_stellar_dex_returns_directly():
    comp = {"Si": 25.0, "Mg": 23.0}
    s = Star(stellar_dex=comp)
    assert s.stellar_dex == comp

# def test_wrong_dex_units_raises():
#     with pytest.raises(ValueError):
#         Star(stellar_dex={"SiO2": 45.0})

def test_unrecognized_units_warns():
    with pytest.warns(UserWarning, match="not recognized"):
        Star(stellar_dex={"Si": 25.0, "Zq": 0.5})

def test_name_returns_directly():
    s = Star(name="zombocom")
    assert s.name == "zombocom"

def test_mass_returns_directly():
    s = Star(mass=1000.1)
    assert s.mass == 1000.1


# ============================================================================
# REVERSE PIPELINE: wtpt_oxides → dex
# ============================================================================
def test_star_from_wtpt_oxides():
    """Star initialized from wtpt_oxides should compute all reverse properties."""
    oxides = {
        'SiO2': 42.37, 'TiO2': 0.07, 'Cr2O3': 0.45,
        'Al2O3': 2.73, 'FeO': 24.25, 'MnO': 0.28,
        'MgO': 25.33, 'NiO': 1.48, 'CaO': 1.69, 'Na2O': 1.35,
    }
    s = Star(wtpt_oxides=oxides)
    assert s.wtpt_oxides is not None
    assert s.wtpt_elements is not None
    assert s.total_wt_atoms is not None
    assert s.atoms_ref_solar is not None
    assert s.ax is not None
    assert s.stellar_dex is not None


def test_star_roundtrip_relative_dex():
    """Star(dex).wtpt_oxides → Star(oxides).stellar_dex preserves relative dex."""
    dex = {
        'Si': 0.27, 'Ti': 4.61436e-16, 'Cr': 0.08, 'Al': 0.23,
        'Fe': 0.02, 'Mn': 0.06, 'Mg': 0.21, 'Ca': 0.1,
        'Na': 0.3, 'Ni': 0.04,
    }
    s1 = Star(stellar_dex=dex)
    oxides = s1.wtpt_oxides

    s2 = Star(wtpt_oxides=oxides)
    recovered = s2.stellar_dex

    ref = 'Fe'
    for el in dex:
        if el == ref:
            continue
        original_diff = dex[el] - dex[ref]
        recovered_diff = recovered[el] - recovered[ref]
        assert recovered_diff == pytest.approx(original_diff, abs=1e-6), \
            f"dex[{el}]-dex[{ref}]: expected {original_diff}, got {recovered_diff}"


def test_star_both_inputs_raises():
    """Passing both stellar_dex and wtpt_oxides should raise ValueError."""
    with pytest.raises(ValueError, match="Cannot pass both"):
        Star(stellar_dex={"Fe": 0.0}, wtpt_oxides={"FeO": 50.0})


def test_star_wtpt_oxides_unrecognized_warns():
    """Unrecognized oxide keys should warn."""
    with pytest.warns(UserWarning, match="not recognized"):
        Star(wtpt_oxides={"SiO2": 45.0, "ZqO": 0.5})
