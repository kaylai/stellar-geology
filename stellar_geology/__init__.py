"""
Stellar Geology

A python library for calculating planet geochemistry from stellar compositions.
"""

__version__ = "0.0.1"
__author__ = "Kayla Iacovino"

from .star import Star
from .planet import Planet
from .constants import *

__all__ = [
    # Classes
    'Star',
    'Planet',
    # Functions
    'createPlanet',
    # Constants
    'oxides_to_elements',
    'elements_to_oxides',
]