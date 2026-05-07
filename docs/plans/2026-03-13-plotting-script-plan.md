# Plotting Script + Exploration Notebook Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extract all plotting logic from `Putirka_and_Rarick_calcs.ipynb` into `sg_plotting.py` and create a clean `exploration.ipynb` notebook where each plot is a single function call with tunable parameters.

**Architecture:** `sg_plotting.py` contains data-generation functions (`load_stellar_csv`, `generate_synthetic_grid`) and one-call plot functions. Each plot function accepts a DataFrame + flexible axis specs (strings evaluated via `df.eval()` or callables). The notebook has two sections: real data (any CSV) and synthetic data (tunable grid params), each followed by plot calls.

**Tech Stack:** Python, matplotlib, numpy, pandas, stellar_geology

**Both files at:** `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/`

---

### Task 1: Create `sg_plotting.py` — Private Helpers

**Files:**
- Create: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/sg_plotting.py`

**Step 1: Write the file with imports and all private helpers**

```python
"""
sg_plotting.py — Plotting functions for stellar geology exploration.

Data generation:
    load_stellar_csv()        — CSV of stellar dex → pipeline → DataFrame
    generate_synthetic_grid() — Synthetic mineralogy grid → DataFrame

Ternary plots:
    plot_mineralogy_ternary()  — Ol-Opx-Cpx projected from garnet
    plot_elemental_ternary()   — Any 3 elements (e.g., Fe-Si-Mg) with optional Mg# isolines
    plot_heatmap_ternary()     — Mineralogy ternary colored by any value (tricontourf)
    plot_bsp_oxide_heatmaps()  — 5-panel BSP oxide heatmaps on mineralogy ternary

Scatter plots:
    plot_scatter()             — Flexible x-y scatter with auto-labels
    plot_bsp_recovery()        — Recovered vs input BSP fidelity (5 panels)

Sweep / line plots:
    plot_ol_opx_sweep()        — Ol/(Ol+Opx) sweep with cpx isolines
    plot_cpx_sweep()           — Cpx sweep with ol/opx ratio isolines
    plot_ca_al_isolines()      — Ca/Al vs Ol/Opx with cpx isolines
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation

import stellar_geology as sg
from stellar_geology import constants as const


# ============================================================
# PRIVATE HELPERS
# ============================================================

_MW_CaO = 56.077
_MW_Al2O3 = 101.961

_PHASE_COLORS = {
    'ol': '#4daf4a',
    'opx': '#377eb8',
    'cpx': '#e41a1c',
    'grt': '#984ea3',
}

_PHASE_LEGEND = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#4daf4a',
           markersize=8, label='Ol-dominated'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#377eb8',
           markersize=8, label='Opx-dominated'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e41a1c',
           markersize=8, label='Cpx-dominated'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#984ea3',
           markersize=8, label='Grt-dominated'),
]


def _eval_column(df, expr):
    """Evaluate a column expression on a DataFrame.

    Parameters
    ----------
    df : DataFrame
    expr : str or callable
        If str, evaluated via df.eval().
        If callable, called as expr(df).

    Returns
    -------
    (Series, str or None)
        The computed values and an auto-label (the string itself, or None
        if a callable was passed).
    """
    if callable(expr):
        return expr(df), None
    return df.eval(expr), str(expr)


def _ternary_to_cartesian(top, bottom_left, bottom_right):
    """Convert ternary coordinates to Cartesian (x, y).

    Vertices:  bottom_left = (0, 0),  bottom_right = (1, 0),
               top = (0.5, sqrt(3)/2).

    Values do NOT need to be positive — negative coordinates extend
    beyond the triangle edges.
    """
    x = bottom_right + 0.5 * top
    y = (np.sqrt(3) / 2) * top
    return x, y


def _draw_ternary_frame(ax, top='Ol', left='Opx', right='Cpx',
                        gridlines=True):
    """Draw equilateral triangle boundary, gridlines, and corner labels."""
    verts = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2], [0, 0]])
    ax.plot(verts[:, 0], verts[:, 1], 'k-', lw=2, zorder=10)

    offset = 0.04
    ax.text(0.5, np.sqrt(3) / 2 + offset, top,
            ha='center', va='bottom', fontsize=14, fontweight='bold')
    ax.text(-offset, -offset, left,
            ha='right', va='top', fontsize=14, fontweight='bold')
    ax.text(1 + offset, -offset, right,
            ha='left', va='top', fontsize=14, fontweight='bold')

    if gridlines:
        for f in np.arange(0.1, 1.0, 0.1):
            # Constant top
            x0, y0 = _ternary_to_cartesian(f, 1 - f, 0)
            x1, y1 = _ternary_to_cartesian(f, 0, 1 - f)
            ax.plot([x0, x1], [y0, y1], 'gray', lw=0.5, alpha=0.4, zorder=1)
            # Constant bottom_left
            x0, y0 = _ternary_to_cartesian(0, f, 1 - f)
            x1, y1 = _ternary_to_cartesian(1 - f, f, 0)
            ax.plot([x0, x1], [y0, y1], 'gray', lw=0.5, alpha=0.4, zorder=1)
            # Constant bottom_right
            x0, y0 = _ternary_to_cartesian(0, 1 - f, f)
            x1, y1 = _ternary_to_cartesian(1 - f, 0, f)
            ax.plot([x0, x1], [y0, y1], 'gray', lw=0.5, alpha=0.4, zorder=1)

    ax.set_aspect('equal')
    ax.axis('off')


def _dominant_phase_color(row):
    """Return color string based on which mineral phase dominates."""
    phases = {'ol': row['ol'], 'opx': row['opx'],
              'cpx': row['cpx'], 'grt': row['grt']}
    dom = max(phases, key=phases.get)
    return _PHASE_COLORS[dom]


def _project_from_garnet(df):
    """Project mineralogy onto Ol-Opx-Cpx face (garnet projected out).

    Returns arrays (xs, ys) in Cartesian coordinates and a boolean mask
    of valid rows.
    """
    ol = df['ol'].values
    opx = df['opx'].values
    cpx = df['cpx'].values
    denom = ol + opx + cpx
    mask = denom != 0
    ol_p = np.where(mask, ol / denom, np.nan)
    opx_p = np.where(mask, opx / denom, np.nan)
    cpx_p = np.where(mask, cpx / denom, np.nan)
    xs, ys = _ternary_to_cartesian(ol_p, opx_p, cpx_p)
    return xs, ys, mask
```

**Step 2: Verify syntax**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "import sg_plotting; print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add sg_plotting.py
git commit -m "feat: sg_plotting.py — imports and private helpers"
```

---

### Task 2: Data Generation Functions

**Files:**
- Modify: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/sg_plotting.py`

**Step 1: Add `generate_synthetic_grid()`**

Append after private helpers:

```python
# ============================================================
# DATA GENERATION
# ============================================================

def generate_synthetic_grid(grid_n=15, grt_fixed=0.03, mg_number=0.89,
                            alphas=None):
    """Generate a synthetic mineralogy grid and run the reverse pipeline.

    Creates a 3-component simplex (ol, opx, cpx) with garnet pinned at
    *grt_fixed*. For each composition, runs:
        mineralogy → BSP → bulk planet → Star → molar abundances

    Parameters
    ----------
    grid_n : int
        Resolution — number of divisions along each simplex edge.
        10 → ~120 pts, 15 → ~133 pts, 20 → ~210 pts.
    grt_fixed : float
        Fixed garnet fraction (default 0.03 = 3%).
    mg_number : float
        Mg# for the reverse mineralogy calculation (default 0.89).
    alphas : dict or None
        Core-partitioning fractions. Default {'Fe': 0.494, 'Ni': 0.08}.

    Returns
    -------
    DataFrame
        Columns: ol, opx, cpx, grt, Si, Fe, Mg, Ca, Al,
        BSP_SiO2, BSP_Al2O3, BSP_FeO, BSP_MgO, BSP_CaO,
        BSP_CaAl_molar, ol_opx_ratio, star_name (NaN).
    """
    if alphas is None:
        alphas = {'Fe': 0.494, 'Ni': 0.08}

    budget = 1.0 - grt_fixed
    atom_mass = const.cationMass

    compositions = []
    for ol_i in range(grid_n + 1):
        for opx_i in range(grid_n + 1 - ol_i):
            cpx_i = grid_n - ol_i - opx_i
            ol = (ol_i / grid_n) * budget
            opx = (opx_i / grid_n) * budget
            cpx = (cpx_i / grid_n) * budget
            compositions.append((ol, opx, cpx))

    # Skip pure end-members
    compositions = [(ol, opx, cpx) for ol, opx, cpx in compositions
                    if max(ol, opx, cpx) < budget * 0.99]

    records = []
    for ol, opx, cpx in compositions:
        mineralogy = {
            'olivine': ol, 'orthopyroxene': opx,
            'clinopyroxene': cpx, 'garnet': grt_fixed,
        }
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bsp = sg.calculate_composition_from_mineralogy(
                    mineralogy, mg_number=mg_number)
                planet = sg.Planet(bulk_silicate_planet=bsp, alphas=alphas)
                bp_oxides = planet.bulk_planet
                star = sg.Star(wtpt_oxides=bp_oxides)
                wt_el = star.wtpt_elements
                if wt_el is None:
                    continue

                moles = {el: wt_el.get(el, 0.0) / atom_mass[el]
                         for el in ['Si', 'Fe', 'Mg', 'Ca', 'Al']}

                rec = {'ol': ol, 'opx': opx, 'cpx': cpx, 'grt': grt_fixed}
                rec.update(moles)
                for ox in ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO']:
                    rec[f'BSP_{ox}'] = bsp.get(ox, 0.0)

                rec['ol_opx_ratio'] = ol / opx if opx > 0 else np.nan
                mol_ca = bsp.get('CaO', 0) / _MW_CaO
                mol_al = 2 * bsp.get('Al2O3', 0) / _MW_Al2O3
                rec['BSP_CaAl_molar'] = (mol_ca / mol_al
                                          if mol_al > 0 else np.nan)
                rec['star_name'] = np.nan
                records.append(rec)
        except Exception:
            continue

    df = pd.DataFrame(records)
    print(f"Generated {len(df)} synthetic compositions "
          f"(grid_n={grid_n}, grt={grt_fixed:.0%})")
    return df
```

**Step 2: Add `load_stellar_csv()`**

```python
def load_stellar_csv(csv_path, alphas=None, mg_number=0.89,
                     dex_columns=None, star_id_column='star_id'):
    """Load a CSV of stellar dex values, run full pipeline, return DataFrame.

    Parameters
    ----------
    csv_path : str
        Path to CSV file.
    alphas : dict or None
        Core-partitioning fractions. Default {'Fe': 0.494, 'Ni': 0.08}.
    mg_number : float
        Mg# for reverse mineralogy (default 0.89).
    dex_columns : list of str or None
        Element columns to read as dex. Default:
        ['Si','Ti','Cr','Al','Fe','Mn','Mg','Ca','Na','Ni'].
    star_id_column : str
        Column name containing star identifiers (default 'star_id').

    Returns
    -------
    DataFrame
        Same schema as generate_synthetic_grid(), plus star_name column.
    """
    if alphas is None:
        alphas = {'Fe': 0.494, 'Ni': 0.08}
    if dex_columns is None:
        dex_columns = ['Si', 'Ti', 'Cr', 'Al', 'Fe', 'Mn', 'Mg', 'Ca',
                        'Na', 'Ni']

    raw = pd.read_csv(csv_path)
    stars_unique = raw.drop_duplicates(subset=star_id_column)
    atom_mass = const.cationMass

    records = []
    for _, row in stars_unique.iterrows():
        star_name = str(row[star_id_column]).strip()
        stellar_dex = {}
        for el in dex_columns:
            val = row.get(el)
            if pd.notna(val):
                stellar_dex[el] = float(val)
        if len(stellar_dex) < 3:
            continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                star = sg.Star(stellar_dex=stellar_dex)
                wt_el_star = star.wtpt_elements
                if wt_el_star is None:
                    continue
                planet = sg.Planet(bulk_planet=wt_el_star,
                                   units="wtpt_elements", alphas=alphas)
                bsp = planet.get_composition(which="bulk_silicate_planet")
                if bsp is None:
                    continue
                mineralogy = sg.calculate_mineralogy(silicate_composition=bsp)

                # Molar abundances from Star
                moles = {el: wt_el_star.get(el, 0.0) / atom_mass.get(el, 1.0)
                         for el in ['Si', 'Fe', 'Mg', 'Ca', 'Al']}

                rec = {
                    'ol': mineralogy.get('olivine', 0),
                    'opx': mineralogy.get('orthopyroxene', 0),
                    'cpx': mineralogy.get('clinopyroxene', 0),
                    'grt': mineralogy.get('garnet', 0),
                }
                rec.update(moles)
                for ox in ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO']:
                    rec[f'BSP_{ox}'] = bsp.get(ox, 0.0)

                ol_val = rec['ol']
                opx_val = rec['opx']
                rec['ol_opx_ratio'] = (ol_val / opx_val
                                        if opx_val > 0 else np.nan)
                mol_ca = bsp.get('CaO', 0) / _MW_CaO
                mol_al = 2 * bsp.get('Al2O3', 0) / _MW_Al2O3
                rec['BSP_CaAl_molar'] = (mol_ca / mol_al
                                          if mol_al > 0 else np.nan)
                rec['star_name'] = star_name
                records.append(rec)
        except Exception as e:
            print(f"Skipped {star_name}: {e}")
            continue

    df = pd.DataFrame(records)
    print(f"Loaded {len(df)} / {len(stars_unique)} stars from {csv_path}")
    return df
```

**Step 3: Verify both functions**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "import sg_plotting as sp; df = sp.generate_synthetic_grid(grid_n=5); print(df.columns.tolist()); print(len(df))"`
Expected: column list printed, ~20 rows

**Step 4: Commit**

```bash
git add sg_plotting.py
git commit -m "feat: add generate_synthetic_grid() and load_stellar_csv()"
```

---

### Task 3: Ternary Plot Functions

**Files:**
- Modify: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/sg_plotting.py`

**Step 1: Add `plot_mineralogy_ternary()`**

```python
# ============================================================
# TERNARY PLOTS
# ============================================================

def plot_mineralogy_ternary(df, labels=False, title=None, ax=None,
                            figsize=(8, 7)):
    """Ol-Opx-Cpx ternary projected from garnet.

    Supports negative cpx (extends beyond the Ol-Opx edge).
    Points colored by dominant mineral phase.

    Parameters
    ----------
    df : DataFrame
        Must contain ol, opx, cpx, grt columns.
    labels : bool
        If True and 'star_name' column exists, annotate each point.
    title : str or None
        Custom title.  Default auto-generates.
    ax : matplotlib Axes or None
        If None, creates a new figure.
    figsize : tuple
        Figure size if creating new figure.

    Returns
    -------
    (fig, ax)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    _draw_ternary_frame(ax, top='Ol', left='Opx', right='Cpx')

    xs, ys, mask = _project_from_garnet(df)
    colors = df.apply(_dominant_phase_color, axis=1).values

    ax.scatter(xs[mask], ys[mask], c=colors[mask], s=30, alpha=0.7,
               edgecolors='k', linewidths=0.3, zorder=5)

    if labels and 'star_name' in df.columns:
        for i in range(len(df)):
            if mask[i] and pd.notna(df.iloc[i].get('star_name')):
                ax.annotate(str(df.iloc[i]['star_name']).strip(),
                            (xs[i], ys[i]), fontsize=6,
                            xytext=(4, 4), textcoords='offset points',
                            alpha=0.7)

    ax.set_xlim(-0.35, 1.15)
    ax.set_ylim(-0.15, 1.05)
    ax.set_title(title or "Ol – Opx – Cpx (projected from Grt)",
                 fontsize=13, pad=15)
    return fig, ax
```

**Step 2: Add `plot_elemental_ternary()`**

```python
def plot_elemental_ternary(df, top='Si', left='Fe', right='Mg',
                           mg_isolines=None, title=None, ax=None,
                           figsize=(8, 7)):
    """Molar-abundance ternary for any 3 elements.

    Parameters
    ----------
    df : DataFrame
        Must contain columns named *top*, *left*, *right* (molar abundances).
    top, left, right : str
        Element names for triangle vertices.
    mg_isolines : list of float or None
        If provided, draw Mg# isolines radiating from *top* vertex.
        Only meaningful when top='Si', left='Fe', right='Mg'.
    title : str or None
    ax : matplotlib Axes or None
    figsize : tuple

    Returns
    -------
    (fig, ax)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    _draw_ternary_frame(ax, top=top, left=left, right=right)

    # Mg# isolines
    if mg_isolines is not None:
        x_top, y_top = _ternary_to_cartesian(1, 0, 0)
        for mgn in mg_isolines:
            x_base, y_base = _ternary_to_cartesian(0, 1 - mgn, mgn)
            ax.plot([x_top, x_base], [y_top, y_base], 'k-', lw=0.8,
                    alpha=0.5, zorder=2)
            ax.text(x_base, y_base - 0.03, f'{mgn:.2f}', ha='center',
                    va='top', fontsize=7, color='k', alpha=0.7)
        ax.text(x_top + 0.12, y_top - 0.05, 'Mg#', fontsize=9,
                fontstyle='italic', color='k', alpha=0.6)

    colors = df.apply(_dominant_phase_color, axis=1).values

    t = df[top].values
    l = df[left].values
    r = df[right].values
    total = t + l + r
    valid = total > 0
    xs, ys = _ternary_to_cartesian(
        np.where(valid, t / total, np.nan),
        np.where(valid, l / total, np.nan),
        np.where(valid, r / total, np.nan),
    )

    ax.scatter(xs[valid], ys[valid], c=colors[valid], s=30, alpha=0.7,
               edgecolors='k', linewidths=0.3, zorder=5)
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.0)
    ax.set_title(title or f"{left} – {top} – {right} (molar)",
                 fontsize=13, pad=15)
    return fig, ax
```

**Step 3: Add `plot_heatmap_ternary()`**

```python
def plot_heatmap_ternary(df, value='BSP_CaAl_molar', cmap='RdYlBu_r',
                         title=None, ax=None, figsize=(8, 7),
                         colorbar_label=None):
    """Mineralogy ternary colored by a value using tricontourf.

    Parameters
    ----------
    df : DataFrame
        Must contain ol, opx, cpx columns and the *value* column.
    value : str or callable
        Column name or expression to color by.
    cmap : str
        Matplotlib colormap name.
    title : str or None
    ax : matplotlib Axes or None
    figsize : tuple
    colorbar_label : str or None
        Label for the colorbar. Auto-generated from *value* if str.

    Returns
    -------
    (fig, ax)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    xs, ys, mask = _project_from_garnet(df)
    vals, auto_label = _eval_column(df, value)
    vals = vals.values if hasattr(vals, 'values') else np.asarray(vals)

    # Drop NaNs
    valid = mask & ~np.isnan(xs) & ~np.isnan(ys) & ~np.isnan(vals)
    xs_c, ys_c, vals_c = xs[valid], ys[valid], vals[valid]

    triang = Triangulation(xs_c, ys_c)
    levels = np.linspace(vals_c.min(), np.quantile(vals_c, 0.95), 30)

    tcf = ax.tricontourf(triang, vals_c, levels=levels, cmap=cmap,
                         extend='both', zorder=2)
    ax.tricontour(triang, vals_c, levels=levels[::3], colors='k',
                  linewidths=0.3, alpha=0.3, zorder=3)

    _draw_ternary_frame(ax)

    cb_label = colorbar_label or auto_label or str(value)
    cb = fig.colorbar(tcf, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label(cb_label, fontsize=12)

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 0.95)
    ax.set_title(title or f"Heatmap: {cb_label}", fontsize=13, pad=10)
    return fig, ax
```

**Step 4: Add `plot_bsp_oxide_heatmaps()`**

```python
def plot_bsp_oxide_heatmaps(df, oxides=None, title=None, figsize=(24, 5)):
    """5-panel BSP oxide heatmaps on the mineralogy ternary.

    Parameters
    ----------
    df : DataFrame
        Must contain ol, opx, cpx and BSP_<oxide> columns.
    oxides : list of str or None
        Which oxides to plot. Default: ['SiO2','Al2O3','FeO','MgO','CaO'].
    title : str or None
    figsize : tuple

    Returns
    -------
    (fig, axes)
    """
    if oxides is None:
        oxides = ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO']
    cmaps = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples']

    xs, ys, mask = _project_from_garnet(df)

    fig, axes = plt.subplots(1, len(oxides), figsize=figsize)
    if len(oxides) == 1:
        axes = [axes]

    for i, ox in enumerate(oxides):
        ax = axes[i]
        _draw_ternary_frame(ax)
        col = f'BSP_{ox}'
        vals = df[col].values
        cmap_name = cmaps[i % len(cmaps)]

        sc = ax.scatter(xs[mask], ys[mask], c=vals[mask], cmap=cmap_name,
                        s=14, alpha=0.7, edgecolors='none', zorder=5)
        cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cb.set_label(f'{ox} (wt%)', fontsize=10)
        ax.set_title(ox, fontsize=13, fontweight='bold')
        ax.set_xlim(-0.35, 1.15)
        ax.set_ylim(-0.15, 1.05)

    fig.suptitle(title or 'BSP Oxide Composition Across Mineralogy Space',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    return fig, axes
```

**Step 5: Verify ternary functions**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "
import sg_plotting as sp
import matplotlib
matplotlib.use('Agg')
df = sp.generate_synthetic_grid(grid_n=5)
fig, ax = sp.plot_mineralogy_ternary(df)
fig.savefig('/tmp/test_min_ternary.png')
fig, ax = sp.plot_elemental_ternary(df, top='Si', left='Fe', right='Mg', mg_isolines=[0.7, 0.89])
fig.savefig('/tmp/test_el_ternary.png')
fig, ax = sp.plot_heatmap_ternary(df, value='BSP_CaAl_molar')
fig.savefig('/tmp/test_heatmap.png')
fig, axes = sp.plot_bsp_oxide_heatmaps(df)
fig.savefig('/tmp/test_bsp_heatmaps.png')
print('All ternary plots OK')
"`
Expected: `All ternary plots OK`

**Step 6: Commit**

```bash
git add sg_plotting.py
git commit -m "feat: add ternary plot functions"
```

---

### Task 4: Scatter Plot Functions

**Files:**
- Modify: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/sg_plotting.py`

**Step 1: Add `plot_scatter()`**

```python
# ============================================================
# SCATTER PLOTS
# ============================================================

def plot_scatter(df, x, y, color_by='phase', xlim=None, ylim=None,
                 xlabel=None, ylabel=None, title=None, ax=None,
                 figsize=(8, 6), s=30, alpha=0.7):
    """General x-y scatter plot with flexible axis expressions.

    Parameters
    ----------
    df : DataFrame
    x, y : str or callable
        Column expression or callable(df) → Series.
        Strings are evaluated via df.eval() and auto-label the axis.
    color_by : str
        'phase' colors by dominant mineral phase (default).
        Any other string is treated as a column name for continuous color.
    xlim, ylim : tuple or None
    xlabel, ylabel : str or None
        Override auto-generated labels.
    title : str or None
    ax : matplotlib Axes or None
    figsize : tuple
    s : float
        Marker size.
    alpha : float
        Marker transparency.

    Returns
    -------
    (fig, ax)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    x_vals, x_auto = _eval_column(df, x)
    y_vals, y_auto = _eval_column(df, y)

    if color_by == 'phase':
        colors = df.apply(_dominant_phase_color, axis=1).values
        ax.scatter(x_vals, y_vals, c=colors, s=s, alpha=alpha,
                   edgecolors='k', linewidths=0.3, zorder=5)
        ax.legend(handles=_PHASE_LEGEND, fontsize=9, loc='best')
    else:
        c_vals = df[color_by].values
        sc = ax.scatter(x_vals, y_vals, c=c_vals, s=s, alpha=alpha,
                        edgecolors='k', linewidths=0.3, zorder=5)
        fig.colorbar(sc, ax=ax, shrink=0.7, label=color_by)

    ax.set_xlabel(xlabel or x_auto or 'x', fontsize=13)
    ax.set_ylabel(ylabel or y_auto or 'y', fontsize=13)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    ax.set_title(title or '', fontsize=13)
    plt.tight_layout()
    return fig, ax
```

**Step 2: Add `plot_bsp_recovery()`**

```python
def plot_bsp_recovery(csv_path, alphas=None, mg_number=0.89,
                      oxides=None, dex_columns=None,
                      star_id_column='star_id', figsize=(20, 4)):
    """Recovered vs input BSP fidelity plot (one panel per oxide).

    Runs the full forward pipeline (dex → BSP → mineralogy) then the
    reverse (mineralogy → recovered BSP) and plots recovered vs input
    with a 1:1 reference line.

    Parameters
    ----------
    csv_path : str
        Path to stellar dex CSV.
    alphas : dict or None
        Default {'Fe': 0.494, 'Ni': 0.08}.
    mg_number : float
        Mg# for reverse calc.
    oxides : list of str or None
        Default ['SiO2','MgO','FeO','CaO','Al2O3'].
    dex_columns : list of str or None
    star_id_column : str
    figsize : tuple

    Returns
    -------
    (fig, axes)
    """
    if alphas is None:
        alphas = {'Fe': 0.494, 'Ni': 0.08}
    if oxides is None:
        oxides = ['SiO2', 'MgO', 'FeO', 'CaO', 'Al2O3']
    if dex_columns is None:
        dex_columns = ['Si', 'Ti', 'Cr', 'Al', 'Fe', 'Mn', 'Mg', 'Ca',
                        'Na', 'Ni']

    raw = pd.read_csv(csv_path)
    stars_unique = raw.drop_duplicates(subset=star_id_column)

    input_bsp = {ox: [] for ox in oxides}
    recov_bsp = {ox: [] for ox in oxides}

    for _, row in stars_unique.iterrows():
        stellar_dex = {}
        for el in dex_columns:
            val = row.get(el)
            if pd.notna(val):
                stellar_dex[el] = float(val)
        if len(stellar_dex) < 3:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                star = sg.Star(stellar_dex=stellar_dex)
                wt_el = star.wtpt_elements
                if wt_el is None:
                    continue
                planet = sg.Planet(bulk_planet=wt_el,
                                   units="wtpt_elements", alphas=alphas)
                bsp = planet.get_composition(which="bulk_silicate_planet")
                if bsp is None:
                    continue
                mineralogy = sg.calculate_mineralogy(
                    silicate_composition=bsp)

                bsp_5 = {k: v for k, v in bsp.items() if k in oxides}
                bsp_5_norm = sg.conversions.normalize_composition(
                    composition=bsp_5, normalization="standard",
                    units="wtpt_oxides")

                rec_bsp = sg.calculate_composition_from_mineralogy(
                    mineralogy, mg_number=mg_number)

                for ox in oxides:
                    input_bsp[ox].append(bsp_5_norm.get(ox, 0))
                    recov_bsp[ox].append(rec_bsp.get(ox, 0))
        except Exception:
            continue

    colors_ox = {'SiO2': '#e41a1c', 'MgO': '#377eb8', 'FeO': '#4daf4a',
                 'CaO': '#984ea3', 'Al2O3': '#ff7f00'}
    markers_ox = {'SiO2': 'o', 'MgO': 's', 'FeO': 'D',
                  'CaO': '^', 'Al2O3': 'v'}

    fig, axes = plt.subplots(1, len(oxides), figsize=figsize, sharey=False)
    if len(oxides) == 1:
        axes = [axes]

    for i, ox in enumerate(oxides):
        ax = axes[i]
        xi = np.array(input_bsp[ox])
        yi = np.array(recov_bsp[ox])
        c = colors_ox.get(ox, 'steelblue')
        m = markers_ox.get(ox, 'o')

        ax.scatter(xi, yi, c=c, marker=m, s=50, edgecolors='k',
                   linewidths=0.5, zorder=5)
        all_vals = np.concatenate([xi, yi])
        if len(all_vals) > 0:
            lo, hi = all_vals.min() * 0.95, all_vals.max() * 1.05
            ax.plot([lo, hi], [lo, hi], 'k--', lw=1, alpha=0.5, label='1:1')
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
        ax.set_xlabel(f'Input BSP {ox} (wt%)', fontsize=11)
        ax.set_ylabel(f'Recovered BSP {ox} (wt%)', fontsize=11)
        ax.set_title(ox, fontsize=13, fontweight='bold', color=c)
        ax.set_aspect('equal')
        ax.legend(fontsize=9, loc='upper left')

    fig.suptitle(f'Reverse Pipeline Fidelity  |  Mg# = {mg_number}',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    return fig, axes
```

**Step 3: Verify**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "
import sg_plotting as sp
import matplotlib
matplotlib.use('Agg')
df = sp.generate_synthetic_grid(grid_n=5)
fig, ax = sp.plot_scatter(df, x='BSP_CaAl_molar', y='ol / (ol + opx)')
fig.savefig('/tmp/test_scatter.png')
print('Scatter OK')
"`
Expected: `Scatter OK`

**Step 4: Commit**

```bash
git add sg_plotting.py
git commit -m "feat: add plot_scatter() and plot_bsp_recovery()"
```

---

### Task 5: Sweep / Line Plot Functions

**Files:**
- Modify: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/sg_plotting.py`

**Step 1: Add `plot_ol_opx_sweep()`**

```python
# ============================================================
# SWEEP / LINE PLOTS
# ============================================================

def plot_ol_opx_sweep(mg_number=0.89, grt_fixed=0.03,
                      cpx_levels=None, oxides=None,
                      n_steps=50, flip=False, title=None,
                      figsize=(22, 4.5)):
    """Sweep Ol/(Ol+Opx) ratio with cpx isolines.

    Generates its own data. Each curve is a different cpx fraction.

    Parameters
    ----------
    mg_number : float
    grt_fixed : float
    cpx_levels : list of float or None
        Cpx fractions for isolines. Default [0.0, 0.05, 0.10, 0.15, 0.20, 0.30].
    oxides : list of str or None
        BSP oxides to plot. Default ['SiO2','Al2O3','FeO','MgO','CaO'].
    n_steps : int
        Number of sweep points per curve.
    flip : bool
        If True, BSP oxide on x-axis and Ol/(Ol+Opx) on y-axis.
    title : str or None
    figsize : tuple

    Returns
    -------
    (fig, axes)
    """
    if cpx_levels is None:
        cpx_levels = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]
    if oxides is None:
        oxides = ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO']

    fig, axes = plt.subplots(1, len(oxides), figsize=figsize,
                             sharey=flip)
    if len(oxides) == 1:
        axes = [axes]

    line_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(cpx_levels)))

    for cpx_val, lc in zip(cpx_levels, line_colors):
        ol_opx_total = 1.0 - cpx_val - grt_fixed
        if ol_opx_total <= 0:
            continue
        ol_fracs = np.linspace(0.01, 0.99, n_steps)
        bsp_results = {ox: [] for ox in oxides}
        ol_frac_valid = []

        for f_ol in ol_fracs:
            ol = f_ol * ol_opx_total
            opx = (1 - f_ol) * ol_opx_total
            mineralogy = {
                'olivine': ol, 'orthopyroxene': opx,
                'clinopyroxene': cpx_val, 'garnet': grt_fixed,
            }
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    bsp = sg.calculate_composition_from_mineralogy(
                        mineralogy, mg_number=mg_number)
                    for ox in oxides:
                        bsp_results[ox].append(bsp.get(ox, 0.0))
                    ol_frac_valid.append(f_ol)
            except Exception:
                continue

        for i, ox in enumerate(oxides):
            if flip:
                axes[i].plot(bsp_results[ox], ol_frac_valid, color=lc,
                             lw=1.8, alpha=0.8, label=f'cpx = {cpx_val:.0%}')
            else:
                axes[i].plot(ol_frac_valid, bsp_results[ox], color=lc,
                             lw=1.8, alpha=0.8, label=f'cpx = {cpx_val:.0%}')

    for i, ox in enumerate(oxides):
        if flip:
            axes[i].set_xlabel(f'BSP {ox} (wt%)', fontsize=11)
            axes[i].set_title(ox, fontsize=13, fontweight='bold')
            axes[i].set_ylim(0, 1)
            if i == 0:
                axes[i].set_ylabel('Ol / (Ol + Opx)', fontsize=11)
        else:
            axes[i].set_xlabel('Ol / (Ol + Opx)', fontsize=11)
            axes[i].set_ylabel(f'BSP {ox} (wt%)', fontsize=11)
            axes[i].set_title(ox, fontsize=13, fontweight='bold')
            axes[i].set_xlim(0, 1)

    axes[-1].legend(fontsize=8, loc='best',
                    title=f'Cpx (grt = {grt_fixed:.0%})',
                    title_fontsize=8)

    default_title = (f'Effect of Ol/Opx Ratio on BSP'
                     f'{"  (flipped)" if flip else ""}'
                     f'  |  Grt = {grt_fixed:.0%}, Mg# = {mg_number}')
    fig.suptitle(title or default_title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig, axes
```

**Step 2: Add `plot_cpx_sweep()`**

```python
def plot_cpx_sweep(mg_number=0.89, grt_fixed=0.03,
                   ol_opx_ratios=None, oxides=None,
                   n_steps=50, title=None, figsize=(22, 4.5)):
    """Sweep cpx fraction with Ol/(Ol+Opx) isolines.

    Parameters
    ----------
    mg_number : float
    grt_fixed : float
    ol_opx_ratios : list of float or None
        Ol/(Ol+Opx) fractions for isolines. Default [0.1, 0.3, 0.5, 0.7, 0.9].
    oxides : list of str or None
    n_steps : int
    title : str or None
    figsize : tuple

    Returns
    -------
    (fig, axes)
    """
    if ol_opx_ratios is None:
        ol_opx_ratios = [0.1, 0.3, 0.5, 0.7, 0.9]
    if oxides is None:
        oxides = ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO']

    fig, axes = plt.subplots(1, len(oxides), figsize=figsize, sharey=False)
    if len(oxides) == 1:
        axes = [axes]

    line_colors = plt.cm.magma(np.linspace(0.15, 0.85, len(ol_opx_ratios)))

    for f_ol, lc in zip(ol_opx_ratios, line_colors):
        cpx_vals = np.linspace(0.01, 1.0 - grt_fixed - 0.02, n_steps)
        bsp_results = {ox: [] for ox in oxides}
        cpx_valid = []

        for cpx_val in cpx_vals:
            ol_opx_total = 1.0 - cpx_val - grt_fixed
            if ol_opx_total <= 0:
                continue
            ol = f_ol * ol_opx_total
            opx = (1 - f_ol) * ol_opx_total
            mineralogy = {
                'olivine': ol, 'orthopyroxene': opx,
                'clinopyroxene': cpx_val, 'garnet': grt_fixed,
            }
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    bsp = sg.calculate_composition_from_mineralogy(
                        mineralogy, mg_number=mg_number)
                    for ox in oxides:
                        bsp_results[ox].append(bsp.get(ox, 0.0))
                    cpx_valid.append(cpx_val)
            except Exception:
                continue

        for i, ox in enumerate(oxides):
            axes[i].plot(cpx_valid, bsp_results[ox], color=lc, lw=1.8,
                         alpha=0.8, label=f'Ol/(Ol+Opx) = {f_ol:.1f}')

    for i, ox in enumerate(oxides):
        axes[i].set_xlabel('Cpx fraction', fontsize=11)
        axes[i].set_ylabel(f'BSP {ox} (wt%)', fontsize=11)
        axes[i].set_title(ox, fontsize=13, fontweight='bold')
        axes[i].set_xlim(0, 1.0 - grt_fixed)

    axes[-1].legend(fontsize=8, loc='best',
                    title=f'Ol/Opx (grt = {grt_fixed:.0%})',
                    title_fontsize=8)

    default_title = (f'Effect of Cpx on BSP  |  '
                     f'Grt = {grt_fixed:.0%}, Mg# = {mg_number}')
    fig.suptitle(title or default_title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig, axes
```

**Step 3: Add `plot_ca_al_isolines()`**

```python
def plot_ca_al_isolines(mg_number=0.89, grt_fixed=0.03,
                        cpx_levels=None, x='BSP_CaAl_molar',
                        y='ol / opx', n_steps=60,
                        xlim=None, ylim=None,
                        xlabel=None, ylabel=None,
                        title=None, ax=None, figsize=(8, 6)):
    """Ca/Al vs Ol/Opx with cpx isolines.

    Generates its own sweep data. Each curve sweeps ol/(ol+opx) at a
    fixed cpx level.

    Parameters
    ----------
    mg_number : float
    grt_fixed : float
    cpx_levels : list of float or None
        Default [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40].
    x : str
        What to plot on x-axis. Supports:
        'BSP_CaAl_molar' (default), or any expression on BSP columns.
    y : str
        What to plot on y-axis. Supports:
        'ol / opx' (default), 'ol / (ol + opx)', etc.
    n_steps : int
    xlim, ylim : tuple or None
    xlabel, ylabel : str or None
    title : str or None
    ax : matplotlib Axes or None
    figsize : tuple

    Returns
    -------
    (fig, ax)
    """
    if cpx_levels is None:
        cpx_levels = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    cpx_colors = plt.cm.plasma(np.linspace(0.15, 0.9, len(cpx_levels)))

    for cpx_val, cc in zip(cpx_levels, cpx_colors):
        ol_opx_total = 1.0 - cpx_val - grt_fixed
        if ol_opx_total <= 0:
            continue
        ol_fracs = np.linspace(0.05, 0.95, n_steps)

        sweep_records = []
        for f_ol in ol_fracs:
            ol = f_ol * ol_opx_total
            opx = (1 - f_ol) * ol_opx_total
            if opx == 0:
                continue
            mineralogy = {
                'olivine': ol, 'orthopyroxene': opx,
                'clinopyroxene': cpx_val, 'garnet': grt_fixed,
            }
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    bsp = sg.calculate_composition_from_mineralogy(
                        mineralogy, mg_number=mg_number)
                    mol_ca = bsp.get('CaO', 0) / _MW_CaO
                    mol_al = 2 * bsp.get('Al2O3', 0) / _MW_Al2O3
                    if mol_al > 0:
                        sweep_records.append({
                            'ol': ol, 'opx': opx, 'cpx': cpx_val,
                            'grt': grt_fixed,
                            'BSP_CaAl_molar': mol_ca / mol_al,
                            'BSP_CaO': bsp.get('CaO', 0),
                            'BSP_Al2O3': bsp.get('Al2O3', 0),
                        })
            except Exception:
                continue

        if not sweep_records:
            continue

        sweep_df = pd.DataFrame(sweep_records)
        x_vals, _ = _eval_column(sweep_df, x)
        y_vals, _ = _eval_column(sweep_df, y)
        ax.plot(x_vals, y_vals, color=cc, lw=2, alpha=0.8,
                label=f'cpx = {cpx_val:.0%}')

    ax.set_xlabel(xlabel or x, fontsize=13)
    ax.set_ylabel(ylabel or y, fontsize=13)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    default_title = (f'Isolines at varying Cpx  |  '
                     f'grt = {grt_fixed:.0%}, Mg# = {mg_number}')
    ax.set_title(title or default_title, fontsize=13)
    ax.legend(fontsize=9, loc='best')
    plt.tight_layout()
    return fig, ax
```

**Step 4: Verify**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "
import sg_plotting as sp
import matplotlib
matplotlib.use('Agg')
fig, axes = sp.plot_ol_opx_sweep(mg_number=0.89)
fig.savefig('/tmp/test_ol_opx.png')
fig, axes = sp.plot_cpx_sweep()
fig.savefig('/tmp/test_cpx.png')
fig, ax = sp.plot_ca_al_isolines(y='ol / (ol + opx)')
fig.savefig('/tmp/test_ca_al.png')
print('All sweep plots OK')
"`
Expected: `All sweep plots OK`

**Step 5: Commit**

```bash
git add sg_plotting.py
git commit -m "feat: add sweep plot functions (ol_opx, cpx, ca_al isolines)"
```

---

### Task 6: Create `exploration.ipynb`

**Files:**
- Create: `/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/exploration.ipynb`

**Step 1: Create notebook with all cells**

The notebook has the following structure:

**Cell 1 (markdown):**
```markdown
# Stellar Geology Exploration Notebook

Single-call plots for investigating mineralogy ↔ composition relationships.

**How to use:** Adjust parameters in the code cells, then run the plot call immediately below. Each plot is self-contained.
```

**Cell 2 (code) — Imports:**
```python
import sg_plotting as sp
import matplotlib.pyplot as plt
```

**Cell 3 (markdown):**
```markdown
---
# Real Data

Load stellar dex values from any CSV. The pipeline runs:
`CSV → Star → Planet → BSP → mineralogy → molar abundances`

**Edit the cell below** to point at your CSV and adjust pipeline parameters.
```

**Cell 4 (code) — Real data params + load:**
```python
# --- EDIT THESE ---
csv_path = "/Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter/hypatia_data/hypatia_P-bearing-stars_09032026.csv"
alphas_real = {'Fe': 0.494, 'Ni': 0.08}
mg_number_real = 0.89

df_real = sp.load_stellar_csv(csv_path, alphas=alphas_real, mg_number=mg_number_real)
df_real.head()
```

**Cell 5 (markdown):**
```markdown
### Mineralogy ternary (real data)
```

**Cell 6 (code):**
```python
sp.plot_mineralogy_ternary(df_real, labels=True, title="Hypatia Stars: Mantle Mineralogy")
plt.show()
```

**Cell 7 (markdown):**
```markdown
### Elemental ternary (real data)
```

**Cell 8 (code):**
```python
# --- EDIT THESE ---
top = 'Si'
left = 'Fe'
right = 'Mg'
mg_isolines = [0.5, 0.6, 0.7, 0.8, 0.89, 0.95]

sp.plot_elemental_ternary(df_real, top=top, left=left, right=right,
                          mg_isolines=mg_isolines)
plt.show()
```

**Cell 9 (markdown):**
```markdown
### Scatter plot (real data)
Change `x_expr` and `y_expr` to explore any relationship.
```

**Cell 10 (code):**
```python
# --- EDIT THESE ---
x_expr = 'BSP_CaAl_molar'
y_expr = 'ol / (ol + opx)'
xlim = (0, 5)
ylim = (0, 1)

sp.plot_scatter(df_real, x=x_expr, y=y_expr, xlim=xlim, ylim=ylim)
plt.show()
```

**Cell 11 (markdown):**
```markdown
### BSP recovery fidelity
How well does the reverse pipeline recover the original BSP?
```

**Cell 12 (code):**
```python
# --- EDIT THESE ---
mg_number_recovery = 0.89

sp.plot_bsp_recovery(csv_path, alphas=alphas_real, mg_number=mg_number_recovery)
plt.show()
```

**Cell 13 (markdown):**
```markdown
---
# Synthetic Data

Generate a grid of mineralogies across Ol-Opx-Cpx space (garnet fixed). The reverse pipeline runs on each point:
`mineralogy → BSP → bulk planet → Star → molar abundances`

**Edit the cell below** to change grid resolution and pipeline parameters.
```

**Cell 14 (code) — Synthetic data params + generate:**
```python
# --- EDIT THESE ---
mg_number = 0.89
alphas = {'Fe': 0.494, 'Ni': 0.08}
grid_n = 15       # 10 = fast, 15 = default, 20+ = smooth
grt_fixed = 0.03  # Fixed garnet fraction

df = sp.generate_synthetic_grid(grid_n=grid_n, grt_fixed=grt_fixed,
                                mg_number=mg_number, alphas=alphas)
df.head()
```

**Cell 15 (markdown):**
```markdown
### Mineralogy ternary (synthetic)
```

**Cell 16 (code):**
```python
sp.plot_mineralogy_ternary(df)
plt.show()
```

**Cell 17 (markdown):**
```markdown
### Elemental ternary (synthetic)
```

**Cell 18 (code):**
```python
# --- EDIT THESE ---
top = 'Si'
left = 'Fe'
right = 'Mg'
mg_isolines = [0.5, 0.6, 0.7, 0.8, 0.89, 0.95]

sp.plot_elemental_ternary(df, top=top, left=left, right=right,
                          mg_isolines=mg_isolines)
plt.show()
```

**Cell 19 (markdown):**
```markdown
### BSP oxide heatmaps on mineralogy ternary
```

**Cell 20 (code):**
```python
sp.plot_bsp_oxide_heatmaps(df)
plt.show()
```

**Cell 21 (markdown):**
```markdown
### Heatmap ternary — any value
Color the mineralogy ternary by any column or expression.
```

**Cell 22 (code):**
```python
# --- EDIT THESE ---
value = 'BSP_CaAl_molar'
cmap = 'RdYlBu_r'

sp.plot_heatmap_ternary(df, value=value, cmap=cmap)
plt.show()
```

**Cell 23 (markdown):**
```markdown
### Scatter plot (synthetic)
Flexible x-y scatter. Change expressions to explore any relationship.
```

**Cell 24 (code):**
```python
# --- EDIT THESE ---
x_expr = 'BSP_CaAl_molar'
y_expr = 'ol / opx'          # Try: 'ol / (ol + opx)', 'cpx', etc.
xlim = (0, 3.0)
ylim = (0, 5.0)

sp.plot_scatter(df, x=x_expr, y=y_expr, xlim=xlim, ylim=ylim)
plt.show()
```

**Cell 25 (markdown):**
```markdown
### Cpx vs Ca/Al
```

**Cell 26 (code):**
```python
sp.plot_scatter(df, x='cpx', y='BSP_CaAl_molar',
                xlim=(0, 1), ylim=(0, df['BSP_CaAl_molar'].quantile(0.98)))
plt.show()
```

**Cell 27 (markdown):**
```markdown
---
## Sweep Plots
These generate their own data internally — edit the parameters below each heading.

### Ol/Opx sweep with cpx isolines
```

**Cell 28 (code):**
```python
# --- EDIT THESE ---
mg_number_sweep = 0.89
grt_fixed_sweep = 0.03
cpx_levels = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30]

sp.plot_ol_opx_sweep(mg_number=mg_number_sweep, grt_fixed=grt_fixed_sweep,
                     cpx_levels=cpx_levels)
plt.show()
```

**Cell 29 (markdown):**
```markdown
### Ol/Opx sweep (flipped axes)
```

**Cell 30 (code):**
```python
sp.plot_ol_opx_sweep(mg_number=mg_number_sweep, grt_fixed=grt_fixed_sweep,
                     cpx_levels=cpx_levels, flip=True)
plt.show()
```

**Cell 31 (markdown):**
```markdown
### Cpx sweep with Ol/(Ol+Opx) isolines
```

**Cell 32 (code):**
```python
# --- EDIT THESE ---
ol_opx_ratios = [0.1, 0.3, 0.5, 0.7, 0.9]

sp.plot_cpx_sweep(mg_number=mg_number_sweep, grt_fixed=grt_fixed_sweep,
                  ol_opx_ratios=ol_opx_ratios)
plt.show()
```

**Cell 33 (markdown):**
```markdown
### Ca/Al vs Ol/Opx with cpx isolines
Change `y_expr` to switch between `ol / opx` and `ol / (ol + opx)`.
```

**Cell 34 (code):**
```python
# --- EDIT THESE ---
y_expr = 'ol / opx'           # Try: 'ol / (ol + opx)'
cpx_isoline_levels = [0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40]
xlim = (0, 5.0)
ylim = (0, 5.0)

sp.plot_ca_al_isolines(mg_number=mg_number_sweep, grt_fixed=grt_fixed_sweep,
                       cpx_levels=cpx_isoline_levels,
                       y=y_expr, xlim=xlim, ylim=ylim)
plt.show()
```

**Step 2: Verify notebook opens**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && python -c "import json; nb = json.load(open('exploration.ipynb')); print(f'{len(nb[\"cells\"])} cells'); print('OK')"`
Expected: `34 cells\nOK`

**Step 3: Commit**

```bash
git add exploration.ipynb
git commit -m "feat: add exploration.ipynb — clean plotting notebook"
```

---

### Task 7: Smoke Test — Run Full Notebook

**Step 1: Run all cells headless**

Run: `cd /Users/kiacovino/Library/CloudStorage/Dropbox/Research/_Projects/_TREC/Jupyter && jupyter nbconvert --to notebook --execute exploration.ipynb --output /tmp/exploration_executed.ipynb --ExecutePreprocessor.timeout=300`

Expected: completes without error, all cells execute.

**Step 2: If errors, fix and re-run**

Fix any import paths, column name mismatches, or API issues.

**Step 3: Final commit**

```bash
git add sg_plotting.py exploration.ipynb
git commit -m "fix: address smoke test issues"
```
