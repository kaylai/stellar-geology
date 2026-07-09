import pandas as pd
import pytest
import plotly.graph_objects as go

from stellar_geology.plot.ternary_plot import (ternary_plot, OverlayConfig,
                                               FigureExporter, DEFAULT_LAYOUT)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
DF = pd.DataFrame({
    'Mg': [30.0, 40.0, 50.0],
    'Si': [30.0, 30.0, 25.0],
    'Fe': [40.0, 30.0, 25.0],
})


def test_returns_figure_with_data():
    fig = ternary_plot(DF, 'Mg', 'Si', 'Fe')
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1
    assert list(fig.data[0].a) == list(DF['Mg'])


def test_default_layout_applied():
    fig = ternary_plot(DF, 'Mg', 'Si', 'Fe')
    assert fig.layout.paper_bgcolor == DEFAULT_LAYOUT['paper_bgcolor']
    assert fig.layout.width == DEFAULT_LAYOUT['width']
    assert fig.layout.font.size == DEFAULT_LAYOUT['font']['size']
    assert fig.layout.ternary.sum == 100
    assert fig.layout.ternary.aaxis.title.text == 'Mg'


def test_layout_kwargs_override_defaults():
    fig = ternary_plot(DF, 'Mg', 'Si', 'Fe',
                       title='My Planets', width=700,
                       font_family='Arial')
    assert fig.layout.title.text == 'My Planets'
    assert fig.layout.width == 700
    assert fig.layout.font.family == 'Arial'
    # Deep merge: overriding font family keeps the default size
    assert fig.layout.font.size == DEFAULT_LAYOUT['font']['size']


def test_missing_column_raises():
    with pytest.raises(ValueError, match="not found in df"):
        ternary_plot(DF, 'Mg', 'Si', 'Xx')


def test_overlay_adds_trace():
    overlay = OverlayConfig(df=DF.iloc[:1], name='Earth')
    fig = ternary_plot(DF, 'Mg', 'Si', 'Fe', overlay=overlay)
    assert len(fig.data) == 2
    assert fig.data[1].name == 'Earth'


def test_overlay_invalid_mode_raises():
    with pytest.raises(ValueError, match="single"):
        OverlayConfig(df=DF, mode='bogus')


def test_exporter_queues_figures():
    exporter = FigureExporter("unused.pdf")
    assert len(exporter) == 0
    exporter.register(ternary_plot(DF, 'Mg', 'Si', 'Fe'))
    exporter.register(ternary_plot(DF, 'Mg', 'Si', 'Fe'))
    assert len(exporter) == 2


def test_exporter_empty_save_raises():
    with pytest.raises(RuntimeError, match="No figures"):
        FigureExporter("unused.pdf").save()


def test_exporter_writes_multipage_pdf(tmp_path):
    pytest.importorskip("kaleido")
    fitz = pytest.importorskip("fitz")

    pdf_path = tmp_path / "figures.pdf"
    exporter = FigureExporter(pdf_path)
    exporter.register(ternary_plot(DF, 'Mg', 'Si', 'Fe'))
    exporter.register(ternary_plot(DF, 'Fe', 'Mg', 'Si'))
    exporter.save()

    with fitz.open(pdf_path) as pdf:
        assert pdf.page_count == 2
    assert len(exporter) == 0  # queue cleared after save
