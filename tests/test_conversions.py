import pytest
from stellar_geology import conversions as conv

# TODO !!
# def test_calculate_molar_mgnum():
#     """wt% oxides to molar Mg# = Mg/(Mg+Fe)"""
#     SILICATE_WTPT_OXIDES = {
#         'SiO2':  50.68,
#         'TiO2':  0.15,
#         'Al2O3': 12.27,
#         'FeO':   6.49,
#         'MnO':   0.12,
#         'MgO':   15.03,
#         'CaO':   13.65,
#         'Na2O':  0.62,
#         'K2O' :  0.17,
#         'P2O5':  0.0,
#         'Cr2O3': 1.06,
#         'NiO':   0.0,
#     }
    
#     EXPECTED_MG_NUM = 0.804975434
#     assert conv.calculate_mg_number(SILICATE_WTPT_OXIDES) == (
#         pytest.approx(EXPECTED_MG_NUM, rel=1e-4))

# ============================================================================
# REVERSE PIPELINE TESTS
# ============================================================================

# Known dex composition of star 32768 from Hypatia catalog
STAR_32768_DEX = {
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

# Known forward result for that star
EXPECTED_BULK_PLANET = {
    'SiO2' : 42.37190416490457,
    'TiO2' : 0.07255698934653063,
    'Cr2O3': 0.44576451274180956,
    'Al2O3': 2.727237231195202,
    'FeO'  : 24.250084830232957,
    'MnO'  : 0.28131691827503635,
    'MgO'  : 25.332083338986063,
    'NiO'  : 1.484564845659013,
    'CaO'  : 1.68697066143704,
    'Na2O' : 1.347516507221789,
}


def test_wtpt_elements_from_wtpt_oxides_roundtrip():
    """oxides → elements → oxides should recover the original composition."""
    elements = conv.calculate_wtpt_elements_from_wtpt_oxides(EXPECTED_BULK_PLANET)
    recovered = conv.calculate_wtpt_oxides_from_wtpt_elements(elements)
    assert recovered == pytest.approx(EXPECTED_BULK_PLANET, rel=1e-10)


def test_dex_roundtrip_relative_values():
    """dex → bulk_planet → dex should recover relative dex values.

    The absolute dex values are shifted by a constant offset (due to the
    normalization step), but the DIFFERENCES between any two elements'
    dex values should be preserved exactly.
    """
    # Forward: dex → oxides
    bulk_planet = conv.calculate_bulk_planet_from_dex(STAR_32768_DEX)

    # Reverse: oxides → dex
    recovered_dex = conv.calculate_dex_from_bulk_planet(bulk_planet)

    # Only compare elements that survive the oxide conversion (no C, O, S)
    common_elements = [el for el in STAR_32768_DEX
                       if el in recovered_dex and el not in ('C', 'O', 'S')]

    # Check that relative dex differences are preserved
    ref_el = 'Fe'
    for el in common_elements:
        if el == ref_el:
            continue
        original_diff = STAR_32768_DEX[el] - STAR_32768_DEX[ref_el]
        recovered_diff = recovered_dex[el] - recovered_dex[ref_el]
        assert recovered_diff == pytest.approx(original_diff, abs=1e-10), \
            f"dex[{el}] - dex[{ref_el}]: expected {original_diff}, got {recovered_diff}"


def test_reverse_pipeline_individual_steps():
    """Each reverse step should be the inverse of its forward counterpart."""
    # Forward pipeline
    ax = conv.calculate_ax_from_dex(STAR_32768_DEX)
    ars = conv.calculate_atoms_ref_solar_from_ax(ax)
    twa = conv.calculate_total_wt_atoms_from_atoms_ref_solar(ars)

    # Reverse each step and check
    ars_recovered = conv.calculate_atoms_ref_solar_from_total_wt_atoms(twa)
    assert ars_recovered == pytest.approx(ars, rel=1e-10)

    ax_recovered = conv.calculate_ax_from_atoms_ref_solar(ars)
    assert ax_recovered == pytest.approx(ax, rel=1e-10)

    dex_recovered = conv.calculate_dex_from_ax(ax)
    # Only check non-volatile elements
    for el in dex_recovered:
        if el not in ('C', 'O', 'S'):
            assert dex_recovered[el] == pytest.approx(STAR_32768_DEX[el], rel=1e-10)


def test_dex_from_ax_warns_on_nonpositive():
    """ax values <= 0 should be skipped with a warning."""
    with pytest.warns(UserWarning, match="non-positive"):
        result = conv.calculate_dex_from_ax({"Fe": 1.0, "Si": 0.0})
    assert "Si" not in result
    assert "Fe" in result
