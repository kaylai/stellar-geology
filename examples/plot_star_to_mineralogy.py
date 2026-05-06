"""
Star to Mineralogy Pipeline
============================

The core workflow: create a star from spectroscopic dex measurements, view its composition, and calculate mantle mineralogy from a bulk silicate planet composition.
"""

# %%
# Create a Star
# -------------
# Stellar compositions are specified in dex notation — logarithmic abundances relative to the Sun. Positive values mean enriched, negative values mean depleted.

import matplotlib.pyplot as plt
import stellar_geology as sg

star = sg.Star(
    stellar_dex={"Fe": 0.1, "Mg": -0.05, "Si": -0.05},
    name="Example Star",
)

# View the star's composition as wt% oxides
composition = star.get_composition(units="wtpt_oxides")
print("Star composition (wt% oxides):")
for oxide, pct in composition.items():
    print(f"  {oxide}: {pct:.2f}%")

# %%
# Calculate Mineralogy
# --------------------
# Given a bulk silicate planet composition (wt% oxides), CIPW normative mineralogy predicts the equilibrium mineral assemblage.

bsp = {"SiO2": 45.0, "Al2O3": 4.0, "FeO": 8.0, "MgO": 38.0, "CaO": 3.5}

mineralogy = sg.calculate_mineralogy(bsp)

print("Mantle mineralogy (molar fractions):")
for mineral, fraction in mineralogy.items():
    print(f"  {mineral}: {fraction:.4f}")

# %%
# Normalize for ternary plotting
# ------------------------------
# The ``plot_norm`` function normalizes olivine, clinopyroxene, and orthopyroxene to sum to 1 — ready for a ternary diagram.

ternary = sg.plot_norm(mineralogy)

print("Ternary-normalized (ol + cpx + opx = 1):")
for mineral, fraction in ternary.items():
    print(f"  {mineral}: {fraction:.4f}")

# %%
# Visualize the mineral assemblage
# ---------------------------------
# A pie chart shows the relative proportions of each mineral phase in the predicted mantle.

# Filter out zero-fraction minerals for a cleaner chart
nonzero = {m: f for m, f in mineralogy.items() if f > 0}
minerals = list(nonzero.keys())
fractions = list(nonzero.values())

cmap = plt.colormaps["Set2"]
colors = [cmap(i) for i in range(len(minerals))]

fig, ax = plt.subplots(figsize=(6, 6))
ax.pie(
    fractions,
    labels=minerals,
    autopct="%.1f%%",
    colors=colors,
    startangle=140,
)
ax.set_title("Mantle Mineralogy\n(CIPW normative, molar fractions)", fontweight="bold")
plt.show()
