import pytest
import csv
import tempfile
import pandas as pd
from pandas.testing import assert_frame_equal
from stellar_geology.batchfile import load_stellar

@pytest.fixture
def sample_csv(tmp_path):
    csv_file = tmp_path / "sample.csv"
    df = pd.DataFrame([
        {'Si': 0.27, 'Ti': 4.61436*10**(-16), 'Cr': 0.08, 'Al': 0.23,
         'Fe': 0.02, 'Mn': 0.06, 'Mg': 0.21, 'Ca': 0.1, 'Na': 0.3, 'Ni': 0.04,
         'C':  -0.14, 'O':  -0.06, 'name': "a_named_star"},
        {'Si': 0.15, 'Ti': 3.6133*10**(-16), 'Cr': 0.01, 'Al': 0.55,
         'Fe': 0.03, 'Mn': 0.02, 'Mg': 0.25, 'Ca': 0.09, 'Na': 0.09, 'Ni': 0.05,
         'C':  -0.19, 'O':  -0.04, 'name': "b_named_star"},
    ])
    df.to_csv(csv_file, index=False)
    return csv_file

def test_load_stellar(sample_csv):
    result = load_stellar(sample_csv)
    
    expected = pd.DataFrame([
        {'Si': 0.27, 'Ti': 4.61436*10**(-16), 'Cr': 0.08, 'Al': 0.23,
         'Fe': 0.02, 'Mn': 0.06, 'Mg': 0.21, 'Ca': 0.1, 'Na': 0.3, 'Ni': 0.04,
         'C':  -0.14, 'O':  -0.06, 'name': "a_named_star"},
        {'Si': 0.15, 'Ti': 3.6133*10**(-16), 'Cr': 0.01, 'Al': 0.55,
         'Fe': 0.03, 'Mn': 0.02, 'Mg': 0.25, 'Ca': 0.09, 'Na': 0.09, 'Ni': 0.05,
         'C':  -0.19, 'O':  -0.04, 'name': "b_named_star"},
    ])
    
    pd.testing.assert_frame_equal(result, expected)
