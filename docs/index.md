# stellar_geology

A Python library for calculating planet geochemistry from stellar compositions.

`stellar_geology` takes a star's elemental abundances (in dex notation),
converts them through a geochemical pipeline, and predicts the bulk and
silicate composition of rocky planets. It then calculates CIPW normative
mineralogy for the resulting silicate mantles — olivine, pyroxene, garnet,
and more.

```python
import stellar_geology as sg

star = sg.Star(stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05, "Ca": 0.02, "Al": 0.03})
planet = sg.Planet.from_star(star)
bsp = sg.Planet(
    bulk_planet=planet.bulk_planet, alphas={"Fe": 0.49, "Ni": 0.49}
).get_composition("bulk_silicate_planet", units="wtpt_oxides")
mineralogy = sg.calculate_mineralogy(bsp)
```

```{toctree}
:maxdepth: 2
:caption: User Guide
:hidden:

getting_started
```

```{toctree}
:maxdepth: 2
:caption: Examples
:hidden:

auto_examples/index
```

```{toctree}
:maxdepth: 2
:caption: API Reference
:hidden:

api
```
