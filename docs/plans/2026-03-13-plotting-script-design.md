# Design: Plotting Script + Clean Exploration Notebook

## Overview

Extract all plotting logic from `Putirka_and_Rarick_calcs.ipynb` into a standalone Python script (`sg_plotting.py`). Create a clean Jupyter notebook (`exploration.ipynb`) where each plot is a single function call with tunable parameters directly above it.

Both files live in `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/`.

## Goals

- One function per plot type, callable from a clean notebook
- Axes are flexible: user passes string expressions (e.g., `'ol / (ol + opx)'`) or callables to control what's plotted
- Same DataFrame schema for real and synthetic data, so all plot functions work with either
- Tunable parameters (Mg#, alphas, grid resolution) sit in notebook cells right above the plot call

## File Layout

| File | Location | Purpose |
|---|---|---|
| `sg_plotting.py` | Jupyter folder | All plot functions + data generation helpers |
| `exploration.ipynb` | Jupyter folder | Clean notebook — tweak params, call one function, see plot |

## Data Generation Functions

```python
load_stellar_csv(csv_path, dex_columns, alphas, mg_number) → DataFrame
```
Reads any CSV of stellar dex values, runs the full forward+reverse pipeline (Star → Planet → BSP → mineralogy → molar abundances), returns a DataFrame.

```python
generate_synthetic_grid(grid_n, grt_fixed, mg_number, alphas) → DataFrame
```
Generates a 3-component simplex grid (ol, opx, cpx) with garnet pinned at `grt_fixed`. Runs reverse pipeline on each composition. Returns DataFrame with identical column schema.

### DataFrame Column Schema

Both functions return a DataFrame with these columns:

| Column | Description |
|---|---|
| `ol`, `opx`, `cpx`, `grt` | Mineral fractions |
| `Si`, `Fe`, `Mg`, `Ca`, `Al` | Molar abundances |
| `BSP_SiO2`, `BSP_Al2O3`, `BSP_FeO`, `BSP_MgO`, `BSP_CaO` | BSP oxide wt% |
| `BSP_CaAl_molar` | Ca/Al molar ratio from BSP |
| `ol_opx_ratio` | ol / opx (NaN when opx=0) |
| `star_name` | Star identifier (real data only, NaN for synthetic) |

## Plot Functions

### Ternary Plots

**`plot_mineralogy_ternary(df, labels=True, title=None, ax=None)`**
Ol-Opx-Cpx projected from garnet. Supports negative cpx (extends beyond Ol-Opx edge). Points colored by dominant phase.

**`plot_elemental_ternary(df, top='Si', left='Fe', right='Mg', mg_isolines=None, title=None, ax=None)`**
Any 3 elements as molar ternary. Optional Mg# isolines radiating from Si vertex.

**`plot_heatmap_ternary(df, value='BSP_CaAl_molar', cmap='RdYlBu_r', title=None, ax=None)`**
Mineralogy ternary filled with tricontourf, colored by any value column or expression.

**`plot_bsp_oxide_heatmaps(df, oxides=None, title=None)`**
5-panel figure (one per BSP oxide) on the mineralogy ternary. Each panel is a scatter heatmap.

### Scatter / Binary Plots

**`plot_bsp_recovery(df_input_bsp, df_recovered_bsp, oxides=None, mg_number=None)`**
Recovered vs input BSP for each oxide. 5 panels with 1:1 reference line.

**`plot_scatter(df, x, y, color_by='phase', xlim=None, ylim=None, xlabel=None, ylabel=None, title=None, ax=None)`**
General x-y scatter. `x` and `y` accept strings (evaluated on DataFrame via `df.eval()`) or callables. Auto-generates axis labels from string expressions. `color_by='phase'` colors by dominant mineral phase.

### Line / Sweep Plots

**`plot_sweep(x, y, isolines, isoline_values, grt_fixed=0.03, mg_number=0.89, alphas=None, n_steps=50, flip=False, xlabel=None, ylabel=None, title=None)`**
Generates its own sweep data internally. Sweeps one variable on x-axis with isolines of another. `x`, `y`, `isolines` are strings identifying the sweep structure.

## Flexible Axis Specification

All `x` and `y` parameters accept:

1. **Strings** — evaluated with `df.eval()`. Supports arithmetic on column names:
   - `'ol / opx'`
   - `'ol / (ol + opx)'`
   - `'BSP_CaO / BSP_Al2O3'`

2. **Callables** — receive DataFrame, return Series:
   - `lambda df: np.log10(df['BSP_CaO'] / 56.077)`

When a string is passed, the axis label auto-generates from it. For callables, use `xlabel`/`ylabel` overrides.

### Internal helper

```python
def _eval_column(df, expr):
    if callable(expr):
        return expr(df), None  # (Series, no auto-label)
    return df.eval(expr), expr  # (Series, string for auto-label)
```

## Private Helpers

| Helper | Purpose |
|---|---|
| `_ternary_to_cartesian(top, bl, br)` | Ternary → Cartesian (x, y) for equilateral triangle |
| `_draw_ternary_frame(ax, top, left, right)` | Triangle boundary + gridlines + labels |
| `_dominant_phase_color(row)` | Color by dominant mineral phase |
| `_eval_column(df, expr)` | String or callable → Series + auto-label |
| `_project_from_garnet(df)` | Renormalize ol+opx+cpx → projected fractions |

## Notebook Structure (`exploration.ipynb`)

```
1. Imports + setup
2. --- REAL DATA ---
   2a. CSV path + pipeline params (editable)
   2b. Load data → df_real
   2c. Plot calls on df_real (one per cell)
3. --- SYNTHETIC DATA ---
   3a. Grid params: mg_number, alphas, grid_n, grt_fixed (editable)
   3b. Generate grid → df_synth
   3c. Plot calls on df_synth (one per cell, params right above)
4. Auto PDF generation (optional final cell)
```

Each plot cell looks like:
```python
# --- Tweak these ---
x_expr = 'BSP_CaAl_molar'
y_expr = 'ol / (ol + opx)'

sg_plotting.plot_scatter(df_synth, x=x_expr, y=y_expr, xlim=(0, 5), ylim=(0, 1))
```
