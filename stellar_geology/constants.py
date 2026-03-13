"""Geochemical constants used throughout stellar_geology.

Includes atomic and oxide molar masses, solar abundance references,
and mapping dictionaries between oxide and element representations.
"""

__all__ = ['oxides_to_elements', 'elements_to_oxides']

cationMass = {
    'Si': 28.0855,
    'Ti': 47.867,
    'Cr': 51.9961,
    'Al': 26.98154,
    'Fe': 55.845,
    'Mn': 54.938,
    'Mg': 24.305,
    'Ni': 58.6934,
    'Ca': 40.078,
    'Na': 22.98977,
    'C' :  12.0107,
    'O' :  15.9994,
    'S' :  32.065,
}

oxideMass = { # volatile-free
    'SiO2' : 60.0835,
    'TiO2' : 79.865,
    'Cr2O3': 151.9892,
    'Al2O3': 101.960077,
    'FeO'  : 71.844,
    'MnO'  : 70.937045,
    'MgO'  : 40.304,
    'NiO'  : 74.6924,
    'CaO'  : 56.077,
    'Na2O' : 61.9785386,
}

oxides_to_elements = { # volatile-free
    'SiO2' : 'Si',
    'TiO2' : 'Ti',
    'Cr2O3': 'Cr',
    'Al2O3': 'Al',
    'FeO'  : 'Fe',
    'MnO'  : 'Mn',
    'MgO'  : 'Mg',
    'NiO'  : 'Ni',
    'CaO'  : 'Ca',
    'Na2O' : 'Na',
}

elements_to_oxides = {v: k for k, v in oxides_to_elements.items()}

A_El = {
    'Si': 7.54,
    'Ti': 4.92,
    'Cr': 5.65,
    'Al': 6.46,
    'Fe': 7.47,
    'Mn': 5.5,
    'Mg': 7.55,
    'Ni': 6.22,
    'Ca': 6.34,
    'Na': 6.3,
    'C' :  8.39,
    'O' :  8.73,
    'S' :  7.16,
}

CationNum = {'SiO2': 1, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 2, 'Na2O': 2,
             'K2O': 2, 'MnO': 1, 'TiO2': 1, 'P2O5': 2, 'Cr2O3': 2,
             'NiO': 1, 'CoO': 1, 'Fe2O3': 2, 'H2O': 2, 'CO2': 1, 'F2O': 2}

OxygenNum = {'SiO2': 2, 'MgO': 1, 'FeO': 1, 'CaO': 1, 'Al2O3': 3, 'Na2O': 1,
             'K2O': 1, 'MnO': 1, 'TiO2': 2, 'P2O5': 5, 'Cr2O3': 3,
             'NiO': 1, 'CoO': 1, 'Fe2O3': 3, 'H2O': 1, 'CO2': 2, 'F2O': 1}