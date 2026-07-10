"""
Stellar Geology

A python library for calculating planet geochemistry from stellar compositions.
"""

from importlib.metadata import version as _version

from .star import Star
from .planet import Planet
from .constants import *
from .mineralogy import calculate_mineralogy, calculate_composition_from_mineralogy, plot_norm
from .plot import ternary_plot

__version__ = _version("stellar_geology")
__author__ = "Kayla Iacovino"

__all__ = [
    # Classes
    'Star',
    'Planet',
    # 'Mineralogy', # to be implemented!
    # Functions
    'calculate_mineralogy',
    'calculate_composition_from_mineralogy',
    'plot_norm',
    # Constants
    'oxides_to_elements',
    'elements_to_oxides',
    # Utils
    'ternary_plot',
]