import pytest
from stellar_geology.mineralogy import (calculate_mineralogy,
                                       calculate_composition_from_mineralogy,
                                       _calculate_mol_prop_ox_cipw,
                                       _calculate_mol_frac_cipw, plot_norm)
from stellar_geology import constants as const

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
XTAL_OXIDES = ["SiO2", "Al2O3", "FeO", "MgO", "CaO"]

# For Putirka and Rarick 2019
BULK_PLANET_OXIDES = {
    'SiO2':  45.0623584204862,
    'TiO2':  0.1692846252333,
    'Al2O3': 3.4813379535577,
    'FeO':   15.3424822627129,
    'MgO':   31.4994657354165,
    'CaO':   2.6937762031481,
    'Na2O':  1.0216136521680,
    'Cr2O3': 0.5905624970884,
    'NiO':   0.1391186501890,
}

EXPECTED_MOL_PROP_OX = {
    'SiO2': 0.749985577,
    'Al2O3': 0.068287638,
    'FeO': 0.213545595,
    'MgO': 0.781539131,
    'CaO': 0.048036753,
}

EXPECTED_MINERALOGY = {
        "olivine": 0.136654901,
        "clinopyroxene": -0.007254847,
        "orthopyroxene": 0.11138965,
        "garnet": 0.026272795,
    }

def test_calculate_mineralogy_Putirka():
    # BSP of C #5 Ringwood from Putirka and Rarick (2019) supplemental spreadsheet 2.
    assert calculate_mineralogy(silicate_composition=BULK_PLANET_OXIDES) == pytest.approx(
        EXPECTED_MINERALOGY, rel=1e-4)

# For Till et al. 2012 JGR
# Composition of B1204 (Till et al. 2012 JGR) in Christy Till's spreadsheet 
# "Plag_Lherz_Regression_v6_Final(in_paper).xlsx" in tab "Plag Lherz All Data"
MANTLE_COMPOSITION_TILL_PLAG = {
    'SiO2':  51.63954907,
    'TiO2':  1.978678745,
    'Al2O3': 17.29911673,
    'FeO':   8.895240639,
    'MnO':   0.155815443,
    'MgO':   6.236804218,
    'CaO':   8.788374359,
    'Na2O':  3.987103774,
    'K2O' :  1.005126937,
    'Cr2O3': 0.014278216,
}

PRESSURE_PLAG = 0.8 # GPa

EXPECTED_MIN_TILL_PLAG = {
    "quartz": 0.000120083399358276,
    "plagioclase": 0.702643038236381,
    "olivine": 0.161793338435394,
    "clinopyroxene": 0.146425441044163,
}

REYKJANES_TILL_PLAG = {
    'SiO2':  48.84,
    'TiO2':  0.74,
    'Al2O3': 17.60,
    'FeO':  7.58,
    'MgO':  11.49,
    'CaO':  12.03,
    'Na2O': 1.68,
    'K2O' : 0.05,
}

PRESSURE_REYKJANES = 0.68

EXPECTED_MIN_REYKJANES = {
    "quartz": 0.044,
    "plagioclase": 0.578,
    "olivine": 0.226,
    "clinopyroxene": 0.152,
}

# Compositoin of INT-A3 (Schwab, B.E., Johnston, A.D., 2001) in Christy Till's spreadsheet 
# "Sp_Lherz_Regression_Final(in_paper).xlsx" in tab "Sp Lherz All Data"
MANTLE_COMPOSITION_TILL_SP = {
    'SiO2':  50.68,
    'TiO2':  0.15,
    'Al2O3': 12.27,
    'FeO':   6.49,
    'MnO':   0.12,
    'MgO':   15.03,
    'CaO':   13.65,
    'Na2O':  0.62,
    'K2O' :  0.17,
    'P2O5':  0.0,
    'Cr2O3': 1.06,
    'NiO':   0.0,
}

PRESSURE_SP = 1 # GPa

EXPECTED_MIN_TILL_SP = {
    "quartz": 0.0896091708197259,
    "plagioclase": 0.375885611588421,
    "olivine": 0.238692442363366,
    "clinopyroxene": 0.295812775228487,
}


##TODO tests for instantiation of the class



def test_calculate_mineralogy_Till_plag():
    assert calculate_mineralogy(silicate_composition=MANTLE_COMPOSITION_TILL_PLAG,
                                pressure = PRESSURE_PLAG, model="Till2012_plag") == pytest.approx(
                                    EXPECTED_MIN_TILL_PLAG, rel=1e-4)

def test_calculate_mineralogy_Till_spin():
    assert calculate_mineralogy(silicate_composition=MANTLE_COMPOSITION_TILL_SP,
                                pressure = PRESSURE_SP, model="Till2012_sp") == pytest.approx(
                                    EXPECTED_MIN_TILL_SP, rel=1e-4)

# TEST INTERMEDIATE CALCS FOR CALCULATION OF MINERALOGY

def test_calculate_mol_prop_ox_cipw():
    assert _calculate_mol_prop_ox_cipw(BULK_PLANET_OXIDES) == (
        pytest.approx(EXPECTED_MOL_PROP_OX, rel=1e-4))    
    
def test_calculate_mol_frac_cipw():
    expected_mol_frac_cipw = {
        'SiO2': 0.402915931515738,
        'Al2O3': 0.036686275447429,
        'FmO': 0.534590932745088,
        'CaO': 0.025806860291746
    }
    assert _calculate_mol_frac_cipw(EXPECTED_MOL_PROP_OX) == (
        pytest.approx(expected_mol_frac_cipw))

def test_plot_norm():
    expected_normed = {
        "olivine": 0.56752800813470,
        "clinopyroxene": -0.03012939065374,
        "orthopyroxene": 0.46260138251904,
    }
    assert plot_norm(EXPECTED_MINERALOGY) == pytest.approx(
        expected_normed, rel=1e-4)

# ---------------------------------------------------------------------------
# Conflicting / invalid inputs raise errors
# ---------------------------------------------------------------------------
def test_plot_norm_insufficient_keys_raises():
    with pytest.raises(ValueError):
        plot_norm({"olivine": 55.0, "clinopyroxene": -2.0})

# ---------------------------------------------------------------------------
# Reverse pipeline: mineralogy → bulk silicate composition
# ---------------------------------------------------------------------------
def test_round_trip_mineralogy():
    """Forward then reverse with correct Mg# recovers the 5 CIPW oxides."""
    mineralogy = calculate_mineralogy(silicate_composition=BULK_PLANET_OXIDES)

    # Compute the true Mg# from the original composition
    mg_mol = BULK_PLANET_OXIDES["MgO"] / const.oxideMass["MgO"]
    fe_mol = BULK_PLANET_OXIDES["FeO"] / const.oxideMass["FeO"]
    mg_number = mg_mol / (mg_mol + fe_mol)

    recovered = calculate_composition_from_mineralogy(mineralogy, mg_number=mg_number)

    # The forward pipeline only uses 5 oxides, so normalize the original to those
    cipw_oxides = ["SiO2", "Al2O3", "FeO", "MgO", "CaO"]
    original_cipw = {k: v for k, v in BULK_PLANET_OXIDES.items() if k in cipw_oxides}
    original_sum = sum(original_cipw.values())
    original_normed = {k: 100.0 * v / original_sum for k, v in original_cipw.items()}

    assert recovered == pytest.approx(original_normed, rel=1e-6)

def test_reverse_mineralogy_missing_phase():
    with pytest.raises(ValueError, match="Missing required mineral phases"):
        calculate_composition_from_mineralogy({"olivine": 0.5, "garnet": 0.1}, mg_number=0.89)

def test_reverse_mineralogy_bad_mg_number():
    with pytest.raises(ValueError, match="mg_number"):
        calculate_composition_from_mineralogy(EXPECTED_MINERALOGY, mg_number=0.0)
    with pytest.raises(ValueError, match="mg_number"):
        calculate_composition_from_mineralogy(EXPECTED_MINERALOGY, mg_number=1.5)


# ---------------------------------------------------------------------------
# Compositions with NaN or missing non-critical oxides
# ---------------------------------------------------------------------------
def test_mineralogy_with_nan_noncritical_oxides():
    """BSP with NaN in non-critical oxides (TiO2, NiO) should still compute
    valid mineralogy after filter strips them."""
    import math
    bsp_with_nan = {
        'SiO2': 50.82, 'TiO2': float('nan'), 'Al2O3': 3.76,
        'Cr2O3': float('nan'), 'FeO': 8.39, 'MnO': float('nan'),
        'MgO': 32.35, 'CaO': 4.4, 'Na2O': 0.27, 'NiO': float('nan'),
    }
    from stellar_geology.planet import Planet
    p = Planet(bulk_silicate_planet=bsp_with_nan)
    result = calculate_mineralogy(p.bulk_silicate_planet)
    assert all(not math.isnan(v) for v in result.values())

def test_mineralogy_with_only_critical_oxides():
    """BSP with only the 5 critical oxides should work fine."""
    bsp_minimal = {'SiO2': 45.0, 'Al2O3': 3.5, 'FeO': 15.3, 'MgO': 31.5, 'CaO': 2.7}
    result = calculate_mineralogy(bsp_minimal)
    assert all(isinstance(v, float) for v in result.values())

def test_mineralogy_missing_oxide_warns_and_fills_zero():
    """BSP missing Al2O3 should warn, fill with 0, and compute valid mineralogy."""
    import math
    bsp_no_al = {'SiO2': 45.0, 'FeO': 15.3, 'MgO': 31.5, 'CaO': 2.7}
    with pytest.warns(UserWarning, match="missing from silicate composition"):
        result = calculate_mineralogy(bsp_no_al)
    assert all(not math.isnan(v) for v in result.values())
