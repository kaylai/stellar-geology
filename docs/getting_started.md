# Getting Started

## Installation

Install from source:

```bash
git clone https://github.com/kaylai/stellar-geology.git
cd stellar-geology
pip install .
```

## Quick Example

Create a star, derive a planet, and calculate its mineralogy:

```python
import stellar_geology as sg

# Create a Star from dex composition
star = sg.Star(
    stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05, "Ca": 0.02, "Al": 0.03},
    name="Example Star",
)

# View the star's composition as wt% oxides
star.get_composition(units="wtpt_oxides")

# Create a Planet from the star
planet = sg.Planet.from_star(star)

# Get the bulk silicate planet composition
# BSP requires core partitioning ratios (alphas) to separate core from mantle
bsp = sg.Planet(
    bulk_planet=planet.bulk_planet, alphas={"Fe": 0.49, "Ni": 0.49}
).get_composition("bulk_silicate_planet", units="wtpt_oxides")

# Calculate CIPW normative mineralogy
mineralogy = sg.calculate_mineralogy(bsp)
# Returns dict with olivine, clinopyroxene, orthopyroxene, garnet fractions
```

This created a star with slightly elevated iron and slightly depleted
magnesium and silicon relative to the Sun (in dex notation), predicted
what a rocky planet orbiting that star would be made of, and calculated
the CIPW normative mineralogy of its silicate mantle. The `alphas`
parameter controls how elements partition between the metallic core
and silicate mantle during differentiation.

## Next Steps

- [Examples Gallery](auto_examples/index) — worked examples showing common workflows
- [API Reference](api) — full documentation of all classes, functions, and constants
