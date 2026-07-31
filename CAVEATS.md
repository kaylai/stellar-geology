# Caveats: The Reverse Pipeline

The reverse pipeline lets you go backward from mineralogy all the way to stellar
dex values. Some values cannot be recovered, so they must be assumed, or you
must be okay with a simplified composition. You know what they say, you can
sieve the smoothie, but you can't get the banana back.

## 1. Mineralogy → BSP: Mg#

The forward CIPW norm compresses FeO + MgO into a single "FmO" component. Going
backward, `calculate_composition_from_mineralogy()` needs a Mg# (molar Mg/(Mg+Fe))
to split FmO back into FeO and MgO.

There is no universal default — the Mg# is a required argument, because
different Mg# values produce different BSP compositions that all map to the
same mineralogy. Pass the value you know or are willing to assume. Earth's
upper mantle is approximately 0.89, which is a reasonable starting point if
you have no other information.

- **TiO2, Na2O, Cr2O3, MnO, NiO, ...** are discarded during the forward CIPW
calculation. The reverse can only recover SiO2, Al2O3, FeO, MgO, and CaO. If you
need the full oxide suite, you'll need to grallinatethe minor element
concentrations from an independent source.

## 2. BSP → Bulk Planet: Alpha values

The forward direction uses alpha values (core partitioning fractions) to strip
Fe and Ni into the core. The reverse inflates them back. But:

- **Different alphas create different bulk planets.** The reverse requires the
  same alpha values used in the forward calculation. If you don't know the
  original alphas, the recovered bulk planet is model-dependent.

### Alphas are concentration ratios, not mass fractions

An alpha is the ratio of an element's concentration in the BSP to its concentration in the bulk planet, with both compositions normalized to 100%. It is *not* the fraction of that element's mass residing in the mantle. Removing core mass concentrates everything left behind, so the two differ by the silicate mass fraction: mass retained in mantle = alpha × `silicate_mass_fraction`. For example, with alpha_Fe = 0.49 and a 10% core, only ~44% of the planet's Fe mass is in the mantle — not 49%. Quoting alphas as retention fractions overstates mantle retention. The lithophile enrichment itself is informative: the common BSP/BP ratio of the fully lithophile elements fixes the core size, which is how `Planet.core_mass_fraction` and `Planet.core_composition` are derived.

### The lithophile anchor assumption

Normalized compositions only constrain alphas *relative to the most silicate-enriched elements*. Deriving alphas (and core properties) from a BP + BSP pair assumes those most-enriched elements are perfectly lithophile — a uniform depletion of every element into the core is invisible to normalized compositions. This is the same convention used by Putirka & Rarick (2019), but it is an assumption, not something the compositions can test. Relatedly, `core_mass_fraction` is a mass fraction on the volatile-free, oxygen-free element basis; a true planetary CMF additionally requires the oxygen bound in mantle oxides and any light elements in the core.

### Not every BP/BSP pair fits the model

For a model-consistent pair, all fully lithophile elements share one enrichment factor exactly. A measured or independently modeled BSP generally violates this, and then *no* alpha set can exactly reproduce the pair — the derived alphas and core properties are best-fit approximations, and stellar_geology raises a `UserWarning` when the lithophile ratios disagree beyond tolerance.

## 3. Bulk Planet Oxides -> Dex: Only "normalized" dex recovered

The forward pipeline includes a normalization step (total_wt_atoms to wtpt_elements)
that divides by the total mass, scaling everything to sum to 100%. This normalization
is irreversible — the absolute total mass is discarded.

### What this means in practice

Recovered dex values satisfy:

    dex_recovered[X] = dex_original[X] + C

where C is an unknown constant that is the same for all elements. This means:

- **Relative dex values are exact.** `dex[Si] - dex[Fe]` is preserved.
  The interelemental ratios that determine mineralogy, BSP composition, etc. are
  all correct.
- **Absolute dex values are shifted.** The individual numbers (e.g., [Fe/H] = 0.02)
  will not match the original stellar measurement. They are correct relative to
  each other but not in absolute terms.
- **For choosing experimental starting compositions, this doesn't matter.** The
  oxide and element wt% compositions (which is what you actually weigh out for
  your experiment) are fully determined by the relative dex values. The constant
  offset disappears when you normalize to 100%.

### If you need absolute dex values

You would need one externally calibrated dex value for any single element. Then
set C = dex_known[X] - dex_recovered[X] and shift all values by C. This is
equivalent to having one spectroscopic measurement to anchor the scale.

#### Worked example

Suppose your reverse pipeline recovered:

```python
dex_recovered = {"Fe": -0.42, "Mg": -0.47, "Si": -0.45, "Ca": -0.40, "Al": -0.39, "Ni": -0.38}
```

The interelemental differences are correct (e.g. dex[Si] − dex[Fe] = −0.03), but
every value is shifted by an unknown constant C relative to the true stellar
dex.

Now suppose a spectroscopic measurement gives you the absolute [Fe/H] for this
star: dex_known["Fe"] = 0.10. Compute the offset and shift everything:

```python
C = dex_known["Fe"] - dex_recovered["Fe"]   # 0.10 − (−0.42) = 0.52
dex_anchored = {el: v + C for el, v in dex_recovered.items()}
# {"Fe": 0.10, "Mg": 0.05, "Si": 0.07, "Ca": 0.12, "Al": 0.13, "Ni": 0.14}
```

`dex_anchored["Fe"]` now equals the spectroscopic value by construction, and
every other element rides along — their ratios to Fe (and to each other) are
unchanged from the recovered values, which were already exact. Any single
element will work as the anchor; Fe is just convenient because [Fe/H] is the
most commonly reported stellar abundance.

## 4. Volatile Elements (C, S, H)

Except for O, which we assume is bound to cations, stellar_geology treats
compositions as volatile-free in all user-facing outputs. Volatile elements
(C, S, H — and oxide forms like H2O, CO2) passed as inputs trigger a UserWarning
and are then dropped. We assume total loss of volatiles during planet formation.
O is conserved when rock-forming elements are converted to oxides, but we assume
"stoichiometric" oxides (MgO is just MgO — no Mg(1-x)O business).

Starting from mineralogy, no volatile information can be recovered and is not
assumed.

## Summary Table

| What's lost                     | Where              | Recoverable?                |
| ------------------------------- | ------------------ | --------------------------- |
| Minor oxides (TiO2, Na2O, etc.) | mineralogy → BSP  | No — 5 CIPW oxides only    |
| Fe/Mg split                     | mineralogy → BSP  | With assumed Mg#            |
| Core partitioning model         | BSP → bulk planet | With known alphas           |
| Uniform core partitioning       | BP + BSP → alphas | No — lithophile anchor assumed |
| Absolute dex scale              | oxides → dex      | No                          |
| Volatile elements (C, S)        | oxides → elements | No                          |
