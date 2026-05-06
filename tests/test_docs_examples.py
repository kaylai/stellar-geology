"""
Test that every code example in the documentation runs without errors.

These tests mirror the exact code shown to users in docs/index.md,
docs/getting_started.md, and the gallery example scripts. If a test
here fails, the corresponding docs snippet is giving users broken code.
"""

import subprocess
import sys
from pathlib import Path

import stellar_geology as sg


# ---- Narrative doc snippets (index.md, getting_started.md) ----

def test_index_md_snippet():
    """Mirrors the hero code block on the landing page (docs/index.md)."""
    star = sg.Star(
        stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05, "Ca": 0.02, "Al": 0.03}
    )
    planet = sg.Planet.from_star(star)
    bsp = sg.Planet(
        bulk_planet=planet.bulk_planet, alphas={"Fe": 0.49, "Ni": 0.49}
    ).get_composition("bulk_silicate_planet", units="wtpt_oxides")
    mineralogy = sg.calculate_mineralogy(bsp)

    assert isinstance(mineralogy, dict)
    assert all(isinstance(v, float) for v in mineralogy.values())


def test_getting_started_snippet():
    """Mirrors the Quick Example in docs/getting_started.md."""
    star = sg.Star(
        stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05, "Ca": 0.02, "Al": 0.03},
        name="Example Star",
    )

    composition = star.get_composition(units="wtpt_oxides")
    assert isinstance(composition, dict)

    planet = sg.Planet.from_star(star)

    bsp = sg.Planet(
        bulk_planet=planet.bulk_planet, alphas={"Fe": 0.49, "Ni": 0.49}
    ).get_composition("bulk_silicate_planet", units="wtpt_oxides")
    assert bsp is not None

    mineralogy = sg.calculate_mineralogy(bsp)
    assert isinstance(mineralogy, dict)
    assert "olivine" in mineralogy


# ---- Gallery example scripts (examples/plot_*.py) ----

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


def _run_example(script_name):
    """Run a gallery script in a subprocess and assert it exits cleanly."""
    script = EXAMPLES_DIR / script_name
    assert script.exists(), f"Example script not found: {script}"
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        env={**__import__("os").environ, "MPLBACKEND": "Agg"},
    )
    assert result.returncode == 0, (
        f"{script_name} failed:\n{result.stderr}"
    )


def test_gallery_star_to_mineralogy():
    """Run examples/plot_star_to_mineralogy.py end-to-end."""
    _run_example("plot_star_to_mineralogy.py")


def test_gallery_unit_conversions():
    """Run examples/plot_unit_conversions.py end-to-end."""
    _run_example("plot_unit_conversions.py")


# def test_gallery_forward_pipeline():
#     """Run examples/plot_forward_pipeline.py end-to-end."""
#     _run_example("plot_forward_pipeline.py")


# def test_gallery_reverse_pipeline():
#     """Run examples/plot_reverse_pipeline.py end-to-end."""
#     _run_example("plot_reverse_pipeline.py")
