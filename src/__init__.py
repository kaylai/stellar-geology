"""
Stellar Geology

A python library for calculating planet geochemistry from stellar compositions.
"""

__version__ = "0.0.1"
__author__ = "Kayla Iacovino"

from .star import Star
from .constants import *

__all__ = [
    # Classes
    'Star',
    # Constants
    'oxides_to_elements',
    'elements_to_oxides',
]