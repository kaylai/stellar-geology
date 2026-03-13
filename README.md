Stellar Geology
===============

[![Documentation Status](https://readthedocs.org/projects/stellar-geology/badge/?version=latest)](https://stellar-geology.readthedocs.io/en/latest/?badge=latest)

A python library for calculating planet geochemistry from stellar compositions.

`stellar_geology` takes a star’s elemental abundances (in dex notation), converts them through a geochemical pipeline, and predicts the bulk and silicate composition of rocky planets. It then calculates CIPW normative mineralogy for the resulting silicate mantles — olivine, pyroxene, garnet, and more.

```python
import stellar_geology as sg

# make a star
star = sg.Star(stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05, "Ca": 0.02, "Al": 0.03})

# derive a planet
planet = sg.Planet.from_star(star)
bsp = sg.Planet(bulk_planet=planet.bulk_planet, alphas={"Fe": 0.49, "Ni": 0.49})
bsp = bsp.get_composition("bulk_silicate_planet", units="wtpt_oxides")

# calculate its mineralogy
mineralogy = sg.calculate_mineralogy(bsp)
```
