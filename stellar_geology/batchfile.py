import pandas as pd

def load_stellar(filepath):
    """
    Loads a csv or excel sheet with rows of star compositions in stellar dex
    notation. Columns should be as elements (e.g., 'Si').
    
    Parameters
    ----------
    filepath    str
        Path to csv or xlsx file to import, with the extension.
    
    Returns
    -------
    pd.DataFrame
        Pandas dataframe object with same layout as input file.  
    """
    # TODO properly implement csv and xls, data cleaning, checks
    df = pd.read_csv(filepath)
    
    # Convert columns to numeric where possible, leaving non-numeric
    # columns (e.g., star names) untouched.
    for col in df.columns:
        converted = pd.to_numeric(df[col], errors='coerce')
        if converted.notna().any():
            df[col] = converted
    
    return df