import pytest
from stellar_geology.planet import Planet
from stellar_geology import conversions as conv

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
BULK_PLANET_OXIDES = {
    'SiO2':  34.3228620623356,
    'TiO2':  0.1289398301558,
    'Al2O3': 2.6516473296260,
    'FeO':   33.7355431278904,
    'MgO':   23.9923487223086,
    'CaO':   2.0517814044420,
    'Na2O':  0.7781373566196,
    'Cr2O3': 0.4498165616991,
    'NiO':   1.8889236049229,
}

ALPHAS = {"Fe": 0.494, "Ni": 0.08, "Si": 0.98}


# ---------------------------------------------------------------------------
# Empty planet
# ---------------------------------------------------------------------------
def test_empty_planet_returns_none():
    p = Planet()
    assert p.bulk_planet is None
    assert p.bulk_silicate_planet is None
    assert p.stellar_dex is None
    assert p.alphas is None
    assert p.name is None
    assert p.mass is None


# ---------------------------------------------------------------------------
# Direct pass-through of init args
# ---------------------------------------------------------------------------
def test_bulk_planet_returns_directly():
    comp = {"SiO2": 45.0, "MgO": 30.0}
    p = Planet(bulk_planet=comp)
    assert p.bulk_planet == comp

def test_bsp_returns_directly():
    bsp = {"SiO2": 50.0, "MgO": 35.0}
    p = Planet(bulk_silicate_planet=bsp)
    assert p.bulk_silicate_planet == bsp

def test_name_returns_directly():
    p = Planet(name="zombocom")
    assert p.name == "zombocom"

def test_mass_returns_directly():
    p = Planet(mass=1000.2)
    assert p.mass == 1000.2

def test_alphas_returns_directly():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas=ALPHAS)
    assert p.alphas == ALPHAS


# ---------------------------------------------------------------------------
# Property accessibility: bulk_planet only (no alphas)
# ---------------------------------------------------------------------------
def test_bulk_planet_only_bp_accessible():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES)
    assert p.bulk_planet is not None
    assert p.bulk_planet == BULK_PLANET_OXIDES

def test_bulk_planet_only_bsp_is_none():
    """BSP cannot be computed without alphas; property returns None silently."""
    p = Planet(bulk_planet=BULK_PLANET_OXIDES)
    assert p.bulk_silicate_planet is None

@pytest.mark.xfail(reason="_calculate_dex_from_bulk not yet implemented")
def test_bulk_planet_only_alphas_is_none():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES)
    assert p.alphas is None


# ---------------------------------------------------------------------------
# Property accessibility: bulk_planet + alphas
# ---------------------------------------------------------------------------
def test_bp_and_alphas_bp_accessible():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas=ALPHAS)
    assert p.bulk_planet is not None

def test_bp_and_alphas_bsp_accessible():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas=ALPHAS)
    bsp = p.bulk_silicate_planet
    assert bsp is not None
    assert isinstance(bsp, dict)
    assert len(bsp) > 0
    assert sum(bsp.values()) == pytest.approx(100.0, rel=1e-4)

def test_bp_and_alphas_alphas_accessible():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas=ALPHAS)
    assert p.alphas == ALPHAS


# ---------------------------------------------------------------------------
# Property accessibility: bulk_silicate_planet only (no alphas)
# ---------------------------------------------------------------------------
def test_bsp_only_bsp_accessible():
    bsp = {"SiO2": 45.0, "MgO": 31.5, "FeO": 15.3, "Al2O3": 3.5, "CaO": 2.7}
    p = Planet(bulk_silicate_planet=bsp)
    assert p.bulk_silicate_planet == bsp

def test_bsp_only_bp_is_none():
    """Bulk planet cannot be computed without alphas; property returns None silently."""
    bsp = {"SiO2": 45.0, "MgO": 31.5, "FeO": 15.3, "Al2O3": 3.5, "CaO": 2.7}
    p = Planet(bulk_silicate_planet=bsp)
    assert p.bulk_planet is None

def test_bsp_only_alphas_is_none():
    bsp = {"SiO2": 45.0, "MgO": 31.5, "FeO": 15.3, "Al2O3": 3.5, "CaO": 2.7}
    p = Planet(bulk_silicate_planet=bsp)
    assert p.alphas is None


# ---------------------------------------------------------------------------
# Property accessibility: name and mass alongside compositions
# ---------------------------------------------------------------------------
def test_name_and_mass_alongside_bulk():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, name="Kepler-42b", mass=0.8)
    assert p.name == "Kepler-42b"
    assert p.mass == 0.8
    assert p.bulk_planet is not None

# ---------------------------------------------------------------------------
# Property accessibility: pipeline for calculated attr
# ---------------------------------------------------------------------------
def test_access_calculated_attr():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas=ALPHAS)
    assert p.bulk_silicate_planet is not None
    assert p.alphas is not None
    assert p.bulk_planet is not None

# ---------------------------------------------------------------------------
# Init with non-default units converts to wtpt_oxides internally
# ---------------------------------------------------------------------------
def test_bulk_planet_with_wtpt_elements_units():
    """Passing bulk_planet in wtpt_elements should convert and store as
    wtpt_oxides, so the .bulk_planet property returns wt% oxides."""
    bp_elements = conv.convert_composition(BULK_PLANET_OXIDES, 'wtpt_elements')
    p = Planet(bulk_planet=bp_elements, units='wtpt_elements')
    assert p.bulk_planet is not None
    # Should have oxide keys, not element keys
    assert 'SiO2' in p.bulk_planet
    assert 'Si' not in p.bulk_planet
    assert p.bulk_planet == pytest.approx(BULK_PLANET_OXIDES, rel=1e-4)

def test_bsp_with_molfrac_oxides_units():
    """Passing BSP in molfrac_oxides should convert and store as wtpt_oxides."""
    bsp_mol = conv.convert_composition(BULK_PLANET_OXIDES, 'molfrac_oxides')
    p = Planet(bulk_silicate_planet=bsp_mol, units='molfrac_oxides')
    assert p.bulk_silicate_planet is not None
    assert 'SiO2' in p.bulk_silicate_planet
    assert p.bulk_silicate_planet == pytest.approx(BULK_PLANET_OXIDES, rel=1e-4)

def test_bulk_planet_with_units_still_computes_bsp():
    """Passing bulk_planet in non-default units + alphas should still
    produce a valid BSP."""
    bp_elements = conv.convert_composition(BULK_PLANET_OXIDES, 'wtpt_elements')
    p = Planet(bulk_planet=bp_elements, alphas=ALPHAS, units='wtpt_elements')
    assert p.bulk_planet is not None
    bsp = p.bulk_silicate_planet
    assert bsp is not None
    assert sum(bsp.values()) == pytest.approx(100.0, rel=1e-4)


# ---------------------------------------------------------------------------
# Conflicting / invalid inputs raise errors
# ---------------------------------------------------------------------------
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

def test_invalid_units_raises():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45.0}, units='invalid_units')


# ---------------------------------------------------------------------------
# Superfluous keys are stripped from compositional dicts
# ---------------------------------------------------------------------------
def test_superfluous_keys_stripped_from_bulk_planet():
    with pytest.warns(UserWarning, match="not recognized"):
        p = Planet(bulk_planet={"SiO2": 45.0, "MgO": 38.0, "Planet": "Earth"})
    assert "Planet" not in p.bulk_planet
    assert p.bulk_planet == {"SiO2": 45.0, "MgO": 38.0}

def test_superfluous_keys_stripped_from_bsp():
    with pytest.warns(UserWarning, match="not recognized"):
        p = Planet(bulk_silicate_planet={"SiO2": 50.0, "MgO": 35.0, "Source": "model"})
    assert "Source" not in p.bulk_silicate_planet
    assert p.bulk_silicate_planet == {"SiO2": 50.0, "MgO": 35.0}

def test_nan_values_replaced_with_zero_in_bulk_planet():
    p = Planet(bulk_planet={"SiO2": 45.0, "MgO": 38.0, "TiO2": float('nan')})
    assert p.bulk_planet["TiO2"] == 0.0
    assert p.bulk_planet == {"SiO2": 45.0, "MgO": 38.0, "TiO2": 0.0}

def test_nan_values_replaced_with_zero_in_bsp():
    p = Planet(bulk_silicate_planet={"SiO2": 50.0, "FeO": 8.0, "MgO": 32.0,
                                     "Al2O3": 3.8, "CaO": 4.4, "NiO": float('nan')})
    assert p.bulk_silicate_planet["NiO"] == 0.0

def test_superfluous_keys_stripped_from_stellar_dex():
    with pytest.warns(UserWarning, match="not recognized"):
        p = Planet(stellar_dex={"Si": 0.27, "Mg": 0.21, "Name": "HD 32768"})
    assert "Name" not in p.stellar_dex

def test_calculate_silicate_from_bulk_noFe():
    with pytest.raises(ValueError):
        Planet(bulk_planet={"SiO2": 45}, alphas=[]).bulk_silicate_planet

# ---------------------------------------------------------------------------
# set_alphas: mutation, cache invalidation, user-input preservation
# ---------------------------------------------------------------------------
def test_set_alphas_recovers_from_forgotten_alphas():
    """If alphas were not passed at construction, set_alphas should let
    bulk_silicate_planet compute on the next access."""
    p = Planet(bulk_planet=BULK_PLANET_OXIDES)  # no alphas
    p.set_alphas(ALPHAS)
    bsp = p.bulk_silicate_planet
    assert bsp is not None
    assert sum(bsp.values()) == pytest.approx(100.0, rel=1e-4)


def test_set_alphas_invalidates_derived_bsp_cache():
    """Changing alphas after a bulk_silicate_planet has been derived must
    cause the next access to recompute with the new alphas, not return the
    stale cached value."""
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas={"Fe": 0.49, "Ni": 0.49})
    bsp1 = dict(p.bulk_silicate_planet)
    p.set_alphas({"Fe": 0.30, "Ni": 0.49})
    bsp2 = p.bulk_silicate_planet
    assert bsp1["FeO"] != pytest.approx(bsp2["FeO"]), (
        "set_alphas did not invalidate the derived bulk_silicate_planet cache"
    )


def test_set_alphas_invalidates_derived_bp_cache():
    """Symmetric: changing alphas must invalidate a bulk_planet that was
    lazily derived from bulk_silicate_planet + previous alphas."""
    bsp = {"SiO2": 45.0, "MgO": 31.5, "FeO": 15.3, "Al2O3": 3.5, "CaO": 2.7}
    p = Planet(bulk_silicate_planet=bsp, alphas={"Fe": 0.49})
    bp1 = dict(p.bulk_planet)
    p.set_alphas({"Fe": 0.30})
    bp2 = p.bulk_planet
    assert bp1["FeO"] != pytest.approx(bp2["FeO"])


def test_set_alphas_preserves_user_provided_bsp():
    """A bulk_silicate_planet supplied by the user is truth, not a cache.
    set_alphas must not clear it."""
    bsp = {"SiO2": 45.0, "MgO": 31.5, "FeO": 15.3, "Al2O3": 3.5, "CaO": 2.7}
    p = Planet(bulk_silicate_planet=bsp, alphas={"Fe": 0.49})
    p.set_alphas({"Fe": 0.30})
    assert p.bulk_silicate_planet == bsp


def test_set_alphas_preserves_user_provided_bp():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas={"Fe": 0.49})
    _ = p.bulk_silicate_planet  # populate derived cache
    p.set_alphas({"Fe": 0.30})
    assert p.bulk_planet == BULK_PLANET_OXIDES


def test_set_alphas_alpha_sweep_produces_distinct_results():
    """Realistic loop: re-using one Planet across an alpha sweep must
    give different BSPs per alpha value."""
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas={"Fe": 0.49, "Ni": 0.49})
    fe_results = []
    for alpha_fe in (0.30, 0.40, 0.49, 0.55):
        p.set_alphas({"Fe": alpha_fe, "Ni": 0.49})
        fe_results.append(p.bulk_silicate_planet["FeO"])
    assert len(set(round(v, 6) for v in fe_results)) == 4


def test_set_alphas_to_none_clears():
    p = Planet(bulk_planet=BULK_PLANET_OXIDES, alphas={"Fe": 0.49})
    p.set_alphas(None)
    assert p.alphas is None or p._alphas is None  # cleared
    with pytest.raises(ValueError, match="alphas is missing"):
        p.get_composition("bulk_silicate_planet", units="wtpt_oxides")
