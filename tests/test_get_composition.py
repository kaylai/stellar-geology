import pytest
from stellar_geology.planet import Planet
from stellar_geology.star import Star
from stellar_geology import conversions as conv

# A known composition used across conversion tests
sample_wtpt_oxides = {
    'SiO2' : 50.0,
    'Al2O3': 15.0,
    'FeO'  : 10.0,
    'MgO'  : 20.0,
    'CaO'  : 5.0,
}


# ============================================================
# UNIT CONVERSION FUNCTION TESTS (convert_composition dispatcher)
# ============================================================

# --- wtpt_oxides ---
def test_convert_wtpt_oxides_identity():
    """wtpt_oxides should return a copy of the input."""
    result = conv.convert_composition(sample_wtpt_oxides, 'wtpt_oxides')
    assert result == pytest.approx(sample_wtpt_oxides)

# --- wtfrac_oxides ---
def test_convert_wtfrac_oxides_sums_to_one():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_oxides')
    assert sum(result.values()) == pytest.approx(1.0)

def test_convert_wtfrac_oxides_is_wtpt_div_100():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_oxides')
    expected = {k: v / 100.0 for k, v in sample_wtpt_oxides.items()}
    assert result == pytest.approx(expected)

def test_convert_wtfrac_oxides_has_oxide_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_oxides')
    assert set(result.keys()) == set(sample_wtpt_oxides.keys())

# --- wtpt_elements ---
def test_convert_wtpt_elements_sums_to_100():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtpt_elements')
    assert sum(result.values()) == pytest.approx(100.0)

def test_convert_wtpt_elements_has_element_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtpt_elements')
    expected_elements = {'Si', 'Al', 'Fe', 'Mg', 'Ca'}
    assert set(result.keys()) == expected_elements

# --- wtfrac_elements ---
def test_convert_wtfrac_elements_sums_to_one():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_elements')
    assert sum(result.values()) == pytest.approx(1.0)

def test_convert_wtfrac_elements_has_element_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_elements')
    expected_elements = {'Si', 'Al', 'Fe', 'Mg', 'Ca'}
    assert set(result.keys()) == expected_elements

def test_convert_wtfrac_elements_is_wtpt_elements_div_100():
    wtpt_el = conv.convert_composition(sample_wtpt_oxides, 'wtpt_elements')
    wtfrac_el = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_elements')
    expected = {k: v / 100.0 for k, v in wtpt_el.items()}
    assert wtfrac_el == pytest.approx(expected)

# --- molfrac_oxides ---
def test_convert_molfrac_oxides_sums_to_one():
    result = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    assert sum(result.values()) == pytest.approx(1.0)

def test_convert_molfrac_oxides_has_oxide_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    assert set(result.keys()) == set(sample_wtpt_oxides.keys())

# --- molfrac_elements ---
def test_convert_molfrac_elements_sums_to_one():
    result = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    assert sum(result.values()) == pytest.approx(1.0)

def test_convert_molfrac_elements_has_element_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    expected_elements = {'Si', 'Al', 'Fe', 'Mg', 'Ca'}
    assert set(result.keys()) == expected_elements

# --- molfrac_singleO ---
def test_convert_molfrac_singleO_has_element_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'molfrac_singleO')
    for key in result:
        assert key not in sample_wtpt_oxides  # keys should be elements, not oxides

# --- molpt_oxides ---
def test_convert_molpt_oxides_sums_to_100():
    result = conv.convert_composition(sample_wtpt_oxides, 'molpt_oxides')
    assert sum(result.values()) == pytest.approx(100.0)

def test_convert_molpt_oxides_has_oxide_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'molpt_oxides')
    assert set(result.keys()) == set(sample_wtpt_oxides.keys())

# --- molpt_elements ---
def test_convert_molpt_elements_sums_to_100():
    result = conv.convert_composition(sample_wtpt_oxides, 'molpt_elements')
    assert sum(result.values()) == pytest.approx(100.0)

def test_convert_molpt_elements_has_element_keys():
    result = conv.convert_composition(sample_wtpt_oxides, 'molpt_elements')
    expected_elements = {'Si', 'Al', 'Fe', 'Mg', 'Ca'}
    assert set(result.keys()) == expected_elements

# --- Relationship tests ---
def test_molpt_oxides_equals_molfrac_oxides_times_100():
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    molpt = conv.convert_composition(sample_wtpt_oxides, 'molpt_oxides')
    expected = {k: v * 100.0 for k, v in molfrac.items()}
    assert molpt == pytest.approx(expected)

def test_molpt_elements_equals_molfrac_elements_times_100():
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    molpt = conv.convert_composition(sample_wtpt_oxides, 'molpt_elements')
    expected = {k: v * 100.0 for k, v in molfrac.items()}
    assert molpt == pytest.approx(expected)

def test_al2o3_doubles_cations():
    """Al2O3 has CationNum=2, so molfrac_elements should have more Al than
    molfrac_oxides has Al2O3, relative to other species."""
    comp = {'SiO2': 50.0, 'Al2O3': 50.0}
    molfrac_ox = conv.convert_composition(comp, 'molfrac_oxides')
    molfrac_el = conv.convert_composition(comp, 'molfrac_elements')
    assert molfrac_el['Al'] / molfrac_el['Si'] > molfrac_ox['Al2O3'] / molfrac_ox['SiO2']

# --- Roundtrip tests ---
def test_molfrac_oxides_roundtrip():
    """Convert wtpt -> molfrac_oxides -> wtpt. Should recover original."""
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    recovered = conv.mol_oxides_to_wtpt_oxides(molfrac)
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_molfrac_elements_roundtrip():
    """Convert wtpt -> molfrac_elements -> wtpt. Should recover original."""
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    recovered = conv.mol_cations_to_wtpt_oxides(molfrac)
    assert recovered == pytest.approx(sample_wtpt_oxides)

# --- Invalid units ---
def test_convert_composition_invalid_units_raises():
    with pytest.raises(ValueError):
        conv.convert_composition(sample_wtpt_oxides, 'invalid_units')


# ============================================================
# NORMALIZATION TESTS
# ============================================================

def test_normalize_standard_wtpt_sums_to_100():
    unnormalized = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0}  # sums to 90
    normed = conv.normalize_composition(unnormalized, normalization='standard',
                                         units='wtpt_oxides')
    assert sum(normed.values()) == pytest.approx(100.0)

def test_normalize_standard_preserves_ratios():
    unnormalized = {'SiO2': 40.0, 'MgO': 20.0}  # 2:1 ratio
    normed = conv.normalize_composition(unnormalized, normalization='standard',
                                         units='wtpt_oxides')
    assert normed['SiO2'] / normed['MgO'] == pytest.approx(2.0)

def test_normalize_standard_mol_sums_to_one():
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    # scale down to make unnormalized
    unnormed = {k: v * 0.5 for k, v in molfrac.items()}
    normed = conv.normalize_composition(unnormed, normalization='standard',
                                         units='molfrac_oxides')
    assert sum(normed.values()) == pytest.approx(1.0)

def test_normalize_fixedvolatiles_preserves_h2o():
    comp = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0, 'H2O': 3.0}
    normed = conv.normalize_composition(comp, normalization='fixedvolatiles',
                                         units='wtpt_oxides')
    assert normed['H2O'] == pytest.approx(3.0)
    non_vol_sum = sum(v for k, v in normed.items() if k not in ('H2O', 'CO2'))
    assert non_vol_sum == pytest.approx(97.0)
    assert sum(normed.values()) == pytest.approx(100.0)

def test_normalize_fixedvolatiles_preserves_h2o_and_co2():
    comp = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0, 'H2O': 3.0, 'CO2': 2.0}
    normed = conv.normalize_composition(comp, normalization='fixedvolatiles',
                                         units='wtpt_oxides')
    assert normed['H2O'] == pytest.approx(3.0)
    assert normed['CO2'] == pytest.approx(2.0)
    non_vol_sum = sum(v for k, v in normed.items() if k not in ('H2O', 'CO2'))
    assert non_vol_sum == pytest.approx(95.0)
    assert sum(normed.values()) == pytest.approx(100.0)

def test_normalize_additionalvolatiles_nonvol_sums_to_100():
    comp = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0, 'H2O': 3.0}
    normed = conv.normalize_composition(comp, normalization='additionalvolatiles',
                                         units='wtpt_oxides')
    non_vol_sum = sum(v for k, v in normed.items() if k not in ('H2O', 'CO2'))
    assert non_vol_sum == pytest.approx(100.0)
    assert normed['H2O'] == pytest.approx(3.0)
    assert sum(normed.values()) == pytest.approx(103.0)

def test_normalize_none_is_identity():
    comp = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0}
    normed = conv.normalize_composition(comp, normalization='none',
                                         units='wtpt_oxides')
    assert normed == comp

def test_normalize_invalid_normalization_raises():
    with pytest.raises(ValueError):
        conv.normalize_composition(sample_wtpt_oxides, normalization='invalid',
                                    units='wtpt_oxides')

def test_normalize_invalid_units_raises():
    with pytest.raises(ValueError):
        conv.normalize_composition(sample_wtpt_oxides, normalization='standard',
                                    units='invalid_units')


# ============================================================
# Star.get_composition() tests
# ============================================================

# known dex composition used across Star tests
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
    'C' : -0.14,
    'O' : -0.06,
}

def test_star_get_composition_wtpt_oxides():
    """get_composition(units='wtpt_oxides') should match the wtptOxides property."""
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='wtpt_oxides')
    assert comp == pytest.approx(s.wtpt_oxides)

def test_star_get_composition_wtpt_elements():
    """get_composition(units='wtpt_elements') should match the wtpt_elements property."""
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='wtpt_elements')
    assert comp == pytest.approx(s.wtpt_elements)

def test_star_get_composition_wtfrac_oxides_sums_to_one():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='wtfrac_oxides')
    assert sum(comp.values()) == pytest.approx(1.0)

def test_star_get_composition_wtfrac_elements_sums_to_one():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='wtfrac_elements')
    assert sum(comp.values()) == pytest.approx(1.0)

def test_star_get_composition_molfrac_oxides_sums_to_one():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='molfrac_oxides')
    assert sum(comp.values()) == pytest.approx(1.0)

def test_star_get_composition_molfrac_elements_sums_to_one():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='molfrac_elements')
    assert sum(comp.values()) == pytest.approx(1.0)

def test_star_get_composition_molfrac_singleO_has_element_keys():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='molfrac_singleO')
    for key in comp:
        assert key in list(star_32768_dex.keys()) or key in ('Mn',)

def test_star_get_composition_molpt_oxides_sums_to_100():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='molpt_oxides')
    assert sum(comp.values()) == pytest.approx(100.0)

def test_star_get_composition_molpt_elements_sums_to_100():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='molpt_elements')
    assert sum(comp.values()) == pytest.approx(100.0)

def test_star_get_composition_with_normalization():
    s = Star(stellar_dex=star_32768_dex)
    comp = s.get_composition(units='wtpt_oxides', normalization='standard')
    assert sum(comp.values()) == pytest.approx(100.0)

def test_star_get_composition_invalid_units_raises():
    s = Star(stellar_dex={'Si': 0.27})
    with pytest.raises(ValueError):
        s.get_composition(units='invalid')

def test_star_get_composition_empty_returns_none():
    s = Star()
    assert s.get_composition(units='wtpt_oxides') is None


# ============================================================
# Planet.get_composition() tests
# ============================================================

def test_planet_get_composition_bulk_planet_wtpt():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='wtpt_oxides')
    assert result == pytest.approx(sample_wtpt_oxides)

def test_planet_get_composition_bulk_planet_wtfrac_oxides():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='wtfrac_oxides')
    assert sum(result.values()) == pytest.approx(1.0)
    expected = {k: v / 100.0 for k, v in sample_wtpt_oxides.items()}
    assert result == pytest.approx(expected)

def test_planet_get_composition_bulk_planet_wtpt_elements():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='wtpt_elements')
    assert sum(result.values()) == pytest.approx(100.0)
    for key in result:
        assert key not in sample_wtpt_oxides  # element keys

def test_planet_get_composition_bulk_planet_wtfrac_elements():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='wtfrac_elements')
    assert sum(result.values()) == pytest.approx(1.0)

def test_planet_get_composition_bulk_planet_molfrac_oxides():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='molfrac_oxides')
    assert sum(result.values()) == pytest.approx(1.0)

def test_planet_get_composition_bulk_planet_molfrac_elements():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='molfrac_elements')
    assert sum(result.values()) == pytest.approx(1.0)

def test_planet_get_composition_bulk_planet_molfrac_singleO():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='molfrac_singleO')
    for key in result:
        assert key not in sample_wtpt_oxides

def test_planet_get_composition_molfrac_singleO_ignores_normalization():
    """Normalization should have no effect on molfrac_singleO results."""
    p = Planet(bulk_planet=sample_wtpt_oxides)
    without = p.get_composition(which='bulk_planet', units='molfrac_singleO')
    with_norm = p.get_composition(which='bulk_planet', units='molfrac_singleO',
                                   normalization='standard')
    assert with_norm == pytest.approx(without)

def test_star_get_composition_molfrac_singleO_normalization_warns():
    """Normalization with molfrac_singleO on Star should warn."""
    s = Star(stellar_dex=star_32768_dex)
    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        s.get_composition(units='molfrac_singleO', normalization='standard')
        assert len(w) >= 1
        assert "not supported" in str(w[-1].message).lower()

def test_star_get_composition_wtpt_elements_normalization_works():
    """Normalization with wtpt_elements on Star should normalize."""
    s = Star(stellar_dex=star_32768_dex)
    result = s.get_composition(units='wtpt_elements', normalization='standard')
    assert sum(result.values()) == pytest.approx(100.0)

def test_planet_get_composition_bulk_planet_molpt_oxides():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='molpt_oxides')
    assert sum(result.values()) == pytest.approx(100.0)

def test_planet_get_composition_bulk_planet_molpt_elements():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    result = p.get_composition(which='bulk_planet', units='molpt_elements')
    assert sum(result.values()) == pytest.approx(100.0)

def test_planet_get_composition_bsp_wtpt():
    comp = {'SiO2': 45.0, 'Al2O3': 10.0, 'FeO': 8.0, 'MgO': 30.0, 'CaO': 7.0}
    p = Planet(bulk_silicate_planet=comp)
    result = p.get_composition(which='bulk_silicate_planet', units='wtpt_oxides')
    assert result == pytest.approx(comp)

def test_planet_get_composition_with_normalization():
    unnormalized = {'SiO2': 40.0, 'MgO': 30.0, 'FeO': 20.0}  # sums to 90
    p = Planet(bulk_planet=unnormalized)
    result = p.get_composition(which='bulk_planet', units='wtpt_oxides',
                                normalization='standard')
    assert sum(result.values()) == pytest.approx(100.0)

def test_planet_get_composition_invalid_which_raises():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    with pytest.raises(ValueError):
        p.get_composition(which='invalid', units='wtpt_oxides')

def test_planet_get_composition_invalid_units_raises():
    p = Planet(bulk_planet=sample_wtpt_oxides)
    with pytest.raises(ValueError):
        p.get_composition(which='bulk_planet', units='invalid')

def test_planet_get_composition_empty_returns_none():
    p = Planet()
    assert p.get_composition(which='bulk_planet', units='wtpt_oxides') is None


# ============================================================
# convert_to_wtpt_oxides roundtrip tests
# ============================================================

def test_convert_to_wtpt_oxides_wtpt_oxides_identity():
    result = conv.convert_to_wtpt_oxides(sample_wtpt_oxides, 'wtpt_oxides')
    assert result == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_wtfrac_oxides():
    converted = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_oxides')
    recovered = conv.convert_to_wtpt_oxides(converted, 'wtfrac_oxides')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_wtpt_elements():
    converted = conv.convert_composition(sample_wtpt_oxides, 'wtpt_elements')
    recovered = conv.convert_to_wtpt_oxides(converted, 'wtpt_elements')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_wtfrac_elements():
    converted = conv.convert_composition(sample_wtpt_oxides, 'wtfrac_elements')
    recovered = conv.convert_to_wtpt_oxides(converted, 'wtfrac_elements')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_molfrac_oxides():
    converted = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    recovered = conv.convert_to_wtpt_oxides(converted, 'molfrac_oxides')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_molfrac_elements():
    converted = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    recovered = conv.convert_to_wtpt_oxides(converted, 'molfrac_elements')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_molfrac_singleO():
    converted = conv.convert_composition(sample_wtpt_oxides, 'molfrac_singleO')
    recovered = conv.convert_to_wtpt_oxides(converted, 'molfrac_singleO')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_molpt_oxides():
    converted = conv.convert_composition(sample_wtpt_oxides, 'molpt_oxides')
    recovered = conv.convert_to_wtpt_oxides(converted, 'molpt_oxides')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_roundtrip_molpt_elements():
    converted = conv.convert_composition(sample_wtpt_oxides, 'molpt_elements')
    recovered = conv.convert_to_wtpt_oxides(converted, 'molpt_elements')
    assert recovered == pytest.approx(sample_wtpt_oxides)

def test_convert_to_wtpt_oxides_invalid_units_raises():
    with pytest.raises(ValueError):
        conv.convert_to_wtpt_oxides(sample_wtpt_oxides, 'invalid')


# ============================================================
# Planet.__init__ with units param
# ============================================================

def test_planet_init_with_wtpt_elements():
    """Create Planet from wtpt_elements, bulk_planet should return wtpt_oxides."""
    elements = conv.convert_composition(sample_wtpt_oxides, 'wtpt_elements')
    p = Planet(bulk_planet=elements, units='wtpt_elements')
    assert p.bulk_planet == pytest.approx(sample_wtpt_oxides)

def test_planet_init_with_molfrac_oxides():
    """Create Planet from molfrac_oxides, bulk_planet should return wtpt_oxides."""
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_oxides')
    p = Planet(bulk_planet=molfrac, units='molfrac_oxides')
    assert p.bulk_planet == pytest.approx(sample_wtpt_oxides)

def test_planet_init_with_molfrac_elements():
    """Create Planet from molfrac_elements, bulk_planet should return wtpt_oxides."""
    molfrac = conv.convert_composition(sample_wtpt_oxides, 'molfrac_elements')
    p = Planet(bulk_planet=molfrac, units='molfrac_elements')
    assert p.bulk_planet == pytest.approx(sample_wtpt_oxides)

def test_planet_init_with_wtpt_oxides_is_identity():
    """Default units='wtpt_oxides' should store as-is."""
    p = Planet(bulk_planet=sample_wtpt_oxides)
    assert p.bulk_planet == pytest.approx(sample_wtpt_oxides)

def test_planet_init_bsp_with_units():
    """bulk_silicate_planet should also be converted when units is specified."""
    bsp = {'SiO2': 45.0, 'Al2O3': 10.0, 'FeO': 8.0, 'MgO': 30.0, 'CaO': 7.0}
    molfrac = conv.convert_composition(bsp, 'molfrac_oxides')
    p = Planet(bulk_silicate_planet=molfrac, units='molfrac_oxides')
    assert p.bulk_silicate_planet == pytest.approx(bsp)

def test_planet_init_invalid_units_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet=sample_wtpt_oxides, units='invalid')
