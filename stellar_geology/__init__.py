"""
Stellar Geology

A python library for calculating planet geochemistry from stellar compositions.
"""

__version__ = "0.0.1"
__author__ = "Kayla Iacovino"

from .star import Star
from .planet import Planet
from .constants import *
from .mineralogy import calculate_mineralogy, calculate_composition_from_mineralogy, plot_norm

__all__ = [
    # Classes
    'Star',
    'Planet',
    # 'Mineralogy', # to be implemented!
    # Functions
    'createPlanet',
    'calculate_mineralogy',
    'calculate_composition_from_mineralogy',
    'plot_norm',
    # Constants
    'oxides_to_elements',
    'elements_to_oxides',
]