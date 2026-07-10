"""
Create a Planet from a Star
===========================

The core workflow: create a star from spectroscopic dex measurements, calculate a derivative bulk silicate planet composition, and visualize its chemistry.
"""

# %%
# Create a Star
# -------------
# Stellar compositions are specified in dex notation — logarithmic abundances relative to the Sun. Positive values mean enriched, negative values mean depleted.

import stellar_geology as sg
import pandas as pd
import plotly.graph_objects as go

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
# Calculate a dervative bulk silicate planet composition
# ------------------------------------------------------
# To get from a `Star()` to a bulk silicate planet composition, we simply create a `Planet()` object from our star, which then allows us to access its bulk composition and bulk silicate composition.
planet = sg.Planet().from_star(star)

# %%
# Let's look at the bulk planet composition first
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

bulk = planet.get_composition(
    which="bulk_planet",
    units="molfrac_elements"
)
# %%
# Now, the bulk _silicate_ composition
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# To do this, we provide a fractionation factor to partition elements between the core and bulk silicate planet via the argument `alphas`. `alphas` is a dictionary that looks like: {"Fe": 0.3, "Ni": 0.1}. In that case, 30% of the planet's Fe and 10% of its Ni will go into the metallic core. "Fe" is required at a minimum, but additional elements may be passed here.

# first let's set the alpha value for our planet
planet.set_alphas({"Fe": 0.3})

# now we can calculate the planet's bulk silicate composition
bsp = planet.get_composition(
    which="bulk_silicate_planet",
    units="molfrac_elements"
)
bsp

# %%
# Now let's visualize these compositions
# --------------------------------------
# First we use stellar_geology's ternary_plot, which is a wrapper for plotly's go.Scatterternary. Any argu ments that you can pass to go.Scatterternary are accepted by sg.ternary_plot and will be passed along. Try passing `ternary=dict(bgcolor="#333333")` as an argument in the ternary_plot() call.

# %%
# .. raw:: html
#
#    <div style="background-color:#fff3cd; border-left:4px solid #ff9800; padding:10px; border-radius:4px">
#    <b>💡 Note:</b> ``ternary_plot`` creates interactive plots, even right here on the ReadTheDocs web page! Try clicking on one of the markers in the legend.
#    </div>
#    &nbsp;


fig = sg.ternary_plot(
   bsp, 
   a="Si",
   b="Fe",
   c="Mg",
   title="Planetary compositions, molar fraction",
   name="BSP",
   width=650,
   base_marker=dict(size=20, color="#ff5252",
                    line=dict(color="#000000", width=1.5)
                    ),
)
fig

# %%
# Add more to the plot
# ^^^^^^^^^^^^^^^^^^^^
# sg.ternary_plot() returns a plotly figure object, so we can add to the plot as we would normally with plotly's go.Scatterternary. Let's add our bulk planet composition.

fig.add_trace(go.Scatterternary(
    a=[bulk["Si"]],
    b=[bulk["Fe"]],
    c=[bulk["Mg"]],
    mode="markers",
    name="Bulk Planet",
    marker=dict(
        size=20,
        color="#ffcd40",
        line=dict(color="#000000", width=1.5)
        )
    ))
fig
