# Configuration file for the Sphinx documentation builder.

import importlib.metadata
import os
import sys

# Add the project root to sys.path so autodoc can find the package
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "stellar_geology"
copyright = "2025, Kayla Iacovino"
author = "Kayla Iacovino"
try:
    version = release = importlib.metadata.version("stellar_geology")
except importlib.metadata.PackageNotFoundError:
    version = release = "unknown"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "sphinx_gallery.gen_gallery",
    "myst_parser",
    "sphinx_design",
]

# -- Options for sphinx-gallery ----------------------------------------------

sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
    "filename_pattern": r"/plot_",
    "remove_config_comments": True,
}

# -- Options for MyST (Markdown support) -------------------------------------

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for autodoc -----------------------------------------------------

autodoc_member_order = "bysource"

# -- Options for Napoleon (NumPy-style docstrings) ---------------------------

napoleon_google_docstrings = False
napoleon_numpy_docstrings = True
napoleon_use_param = False
napoleon_use_ivar = True

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_logo = "img/logo.png"
html_theme_options = {
    "logo_only": True,
    "version_selector": True,
    "style_nav_header_background": "#1a1147",
}

# Syntax highlighting style
pygments_style = "manni"

# -- Intersphinx mapping -----------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- Nitpick mode (catch broken references) ----------------------------------

nitpicky = True
nitpick_ignore = [
    # Add entries here as needed, e.g.:
    # ("py:class", "numpy.ndarray"),
]

