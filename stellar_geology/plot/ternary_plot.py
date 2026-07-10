"""
ternary_plot.py
---------------
Lightweight wrapper around Plotly's ternary scatter for consistent
styling and optional data overlays.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_BASE_MARKER = dict(
    size=6, color="black", symbol="circle",
    line=dict(color="black", width=0.5)
)

DEFAULT_OVERLAY_MARKER = dict(
    size=8, color="red", symbol="circle"
)

DEFAULT_AXIS_STYLE  = dict(
    color="black",
    gridcolor="black",
    linecolor="black",
    linewidth=1,
    )

DEFAULT_LAYOUT = dict(
    paper_bgcolor="rgba(0, 0, 0, 0)",
    font=dict(size=18),
    width=900,
    height=800,
    )


# ---------------------------------------------------------------------------
# Overlay configuration
# ---------------------------------------------------------------------------

@dataclass
class OverlayConfig:
    """
    Configuration for overlaying additional data on a ternary plot.

    Parameters
    ----------
    df : pd.DataFrame
        Data to overlay.
    mode : {"single", "per_row"}
        ``"single"`` plots all rows as one trace; ``"per_row"`` gives each
        row its own trace and legend entry (row index label).
    name : str
        Legend label. Used in ``"single"`` mode only.
    marker : dict or list of dict or None
        Marker style. Single dict for ``"single"`` mode; list of dicts
        (cycled) for ``"per_row"`` mode. If None, a default red circle
        is used.
    """
    df: pd.DataFrame
    mode: str = "single"
    name: str = "Overlay"
    marker: dict | list[dict] | None = None

    def __post_init__(self) -> None:
        if self.mode not in ("single", "per_row"):
            raise ValueError(
                f"OverlayConfig.mode must be 'single' or 'per_row', got '{self.mode}'"
            )


# ---------------------------------------------------------------------------
# Figure exporter
# ---------------------------------------------------------------------------

class FigureExporter:
    """
    Collect Plotly figures and write them to a multi-page PDF.

    Parameters
    ----------
    output_path : str or os.PathLike
        Destination path for the output PDF file.

    Examples
    --------
    >>> exporter = FigureExporter("figures.pdf")
    >>> fig = ternary_plot(...)
    >>> exporter.register(fig)
    >>> exporter.save()
    """

    def __init__(self, output_path: str | os.PathLike) -> None:
        self.output_path = output_path
        self._queue: list[go.Figure] = []

    def register(self, fig: go.Figure) -> None:
        """
        Queue a figure for PDF export.

        Parameters
        ----------
        fig : go.Figure
            The figure to queue.
        """
        self._queue.append(fig)
        print(f"Figure queued ({len(self._queue)} total). Call save() to write PDF.")

    def save(self) -> None:
        """
        Write all queued figures to a multi-page PDF and clear the queue.

        Raises
        ------
        RuntimeError
            If the queue is empty.
        """
        import fitz
        import plotly.io as pio

        if not self._queue:
            raise RuntimeError("No figures in queue. Call register() first.")

        merged_pdf = fitz.open()
        for fig in self._queue:
            img_bytes = pio.to_image(fig, format="pdf")
            single_pdf = fitz.open("pdf", img_bytes)
            merged_pdf.insert_pdf(single_pdf)

        merged_pdf.save(self.output_path)
        merged_pdf.close()
        print(f"PDF saved → {self.output_path} ({len(self._queue)} figures)")
        self._queue.clear()

    def __len__(self) -> int:
        return len(self._queue)

    def __repr__(self) -> str:
        return f"FigureExporter(output_path={self.output_path!r}, queued={len(self._queue)})"


# ---------------------------------------------------------------------------
# Ternary plot
# ---------------------------------------------------------------------------

def ternary_plot(
        df: pd.DataFrame | dict,
        a: str,
        b: str,
        c: str,
        name: str | None = None,
        base_marker: dict | None = None,
        axis_style: dict | None = None,
        overlay: OverlayConfig | None = None,
        aaxis_min: float | None = None,
        **layout_kwargs: Any,
) -> go.Figure:
    """
    Create a ternary scatter plot, optionally overlaying additional data.

    Parameters
    ----------
    df : pd.DataFrame or dict
        Primary dataset to plot. A dict may map column names to lists
        (multiple points) or to scalars (a single point), e.g.
        ``{'Mg': 30, 'Si': 40, 'Fe': 30}``.
    a : str
        Column name for the top axis component.
    b : str
        Column name for the left axis component.
    c : str
        Column name for the right axis component.
    name : str, optional
        Legend label for the primary dataset trace. Default is ``"Data"``.
    base_marker : dict, optional
        Plotly marker dict for the primary trace. If None, defaults to a
        small black circle with a black outline.
    axis_style : dict, optional
        Plotly axis formatting kwargs unpacked into each ternary axis (e.g.
        tickfont, gridcolor). If None, default styling is applied.
    overlay : OverlayConfig, optional
        Overlay configuration. If None, no overlay is drawn.
    aaxis_min : float, optional
        Minimum value for the a-axis. If None, Plotly uses its default.
    **layout_kwargs
        Any Plotly layout property, forwarded to ``fig.update_layout()``
        after the defaults in ``DEFAULT_LAYOUT`` are applied (transparent
        background, font size 18, 900x800 figure), so user values win.
        See https://plotly.com/python/reference/layout/ for all options.
        Plotly's magic-underscore shorthand works (e.g. ``font_family``).

    Returns
    -------
    fig : go.Figure
        The constructed ternary scatter figure.

    Examples
    --------
    Commonly used layout kwargs:

    >>> fig = ternary_plot(
    ...     df, "Mg", "Si", "Fe",
    ...     title="Planet compositions",
    ...     width=700,
    ...     height=600,
    ...     paper_bgcolor="white",
    ...     font_size=14,
    ...     showlegend=False,
    ...     legend_title_text="Dataset",
    ... )
    
    You can also pass fig.update_layout:
    
    >>> fig.update_layout(ternary=dict(bgcolor="#676767"))

    Raises
    ------
    ValueError
        If `a`, `b`, or `c` are not columns in `df`.
    """
    if isinstance(df, dict):
        try:
            df = pd.DataFrame(df)
        except ValueError:
            # dict of all scalars (a single point) needs an explicit index
            df = pd.DataFrame(df, index=[0])

    for col in (a, b, c):
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' not found in df. "
                f"Available columns: {df.columns.tolist()}"
            )

    if base_marker is None:
        base_marker = DEFAULT_BASE_MARKER
    if axis_style is None:
        axis_style = DEFAULT_AXIS_STYLE

    fig = px.scatter_ternary(df, a=a, b=b, c=c)
    fig.update_traces(
        marker=base_marker,
        name=name,
        showlegend=True,
        selector=dict(mode="markers"),
    )

    if overlay is not None:
        _add_overlay(fig, overlay, a, b, c)

    axis_kwargs = {}
    if aaxis_min is not None:
        axis_kwargs["aaxis_min"] = aaxis_min

    fig.update_layout(
        ternary=dict(
            sum=100,
            aaxis=dict(title=a, **axis_style),
            baxis=dict(title=b, **axis_style),
            caxis=dict(title=c, **axis_style),
            bgcolor="#F8F8F8",
            **axis_kwargs,
        ),
        **DEFAULT_LAYOUT,
    )
    # Applied second so user kwargs deep-merge over the defaults
    fig.update_layout(**layout_kwargs)
    return fig


def _add_overlay(
        fig: go.Figure,
        overlay: OverlayConfig,
        a: str,
        b: str,
        c: str,
) -> None:
    """Add overlay traces to an existing ternary figure."""
    if overlay.mode == "single":
        marker = overlay.marker if isinstance(overlay.marker, dict) else DEFAULT_OVERLAY_MARKER
        fig.add_trace(go.Scatterternary(
            a=overlay.df[a].tolist(),
            b=overlay.df[b].tolist(),
            c=overlay.df[c].tolist(),
            mode="markers",
            name=overlay.name,
            marker=marker,
        ))

    else:  # per_row
        markers = overlay.marker if isinstance(overlay.marker, list) else [DEFAULT_OVERLAY_MARKER]
        for i, (idx, row) in enumerate(overlay.df.iterrows()):
            fig.add_trace(go.Scatterternary(
                a=[row[a]], b=[row[b]], c=[row[c]],
                mode="markers",
                name=str(idx),
                marker=markers[i % len(markers)],
            ))
