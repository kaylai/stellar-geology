# Caveats: The Reverse Pipeline

The reverse pipeline lets you go backward from mineralogy all the way to stellar
dex values. Some values cannot be recovered, so they must be assumed, or you
must be okay with a simplified composition. You know what they say, you can
sieve the smoothie, but you can't get the banana back.

## 1. Mineralogy → BSP: Mg#

The forward CIPW norm compresses FeO + MgO into a single "FmO" component. Going
backward, `calculate_composition_from_mineralogy()` needs a Mg# (molar Mg/(Mg+Fe))
to split FmO back into FeO and MgO.

The default is 0.89, which is approximately Earth's upper mantle.Different Mg#
values will produce different BSP compositions that all map to the same
mineralogy. If you're working with a composition where you know the Mg#, pass it
explicitly.

- **TiO2, Na2O, Cr2O3, MnO, NiO, ...** are discarded during the forward CIPW
calculation. The reverse can only recover SiO2, Al2O3, FeO, MgO, and CaO. If you
need the full oxide suite, you'll need to grallinatethe minor element
concentrations from an independent source.

## 2. BSP → Bulk Planet: Alpha values

The forward direction uses alpha values (core partitioning fractions) to strip
Fe and Ni into the core. The reverse inflates them back. But:

- **Different alphas → different bulk planets.** The reverse requires the
  *same* alpha values used in the forward calculation. If you don't know the
  original alphas, the recovered bulk planet is model-dependent.

## 3. Bulk Planet Oxides → Dex: Only "normalized" dex recovered

The forward pipeline includes a normalization step (total_wt_atoms → wtpt_elements)
that divides by the total mass, scaling everything to sum to 100%. This normalization
is *irreversible* — the absolute total mass is discarded.

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

## 4. Volatile Elements (C, O, S)

The forward pipeline carries volatile elements (C, O, S) through unit conversion
steps (they appear in the element wt% output). We assume total loss of elements
like C, S, and H. O is conserved when elements are converted to oxides, but we
assume "stoichiometric" oxides (MgO is just MgO, though. No 1-x business).

Starting from mineralogy, no volatile information can be recovered and is not
assumed.

## Summary Table

| What's lost                     | Where              | Recoverable?                |
| ------------------------------- | ------------------ | --------------------------- |
| Minor oxides (TiO2, Na2O, etc.) | mineralogy → BSP  | No — 5 CIPW oxides only    |
| Fe/Mg split                     | mineralogy → BSP  | With assumed Mg#            |
| Core partitioning model         | BSP → bulk planet | With known alphas           |
| Absolute dex scale              | oxides → dex      | No                          |
| Volatile elements (C, S)        | oxides → elements | No                          |
