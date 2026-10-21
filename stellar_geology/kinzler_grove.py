"""
Equation to calculate the Mg# defined as Mg/(Mg+Fe) with total Fe as Fe2+ from Kinzler and Grove
(1992) Primary Magmas of Mid-Ocean Ridge Basalts, Journal of Geophysical Research, vol. 97, no. B5,
pp. 6885–6909.

Mg#m = (
    [MgO]b - (1-F) * (
        Xr_ol * (66.67 * Mg#m / (Mg#m + KD_olm * (1-Mg#m))) + 
        Xr_opx * (50.0 * Mg#m / (Mg#m + KD_opxm * (1-Mg#m))) +
        Xr_aug * (25.0) * Mg#m / (Mg#m + KD_augm * (1-Mg#m))) +
        Xr_spi * (33.33 * Mg#m / (Mg#m + KD_spim * (1-Mg#m)))
    ) / (
    [MgO]b + [FeO]b - (1-F) * (
        Xr_ol * 66.67 + Xr_opx * 50.0 + Xr_aug * 25.0 + Xr_spi * 33.33
    ))
)

Mg#m: Mg# of the melt
[i]b: Cation % of component i in initial mantle composition
F: Melt fraction out of 1
Xr_xt: Phase proportion of mineral xt left after melt is produced, mole fraction
KD_xtm: Fe-Mg crystal melt exchange coefficient of crystal xt
ol: olivine
opx: orthopyroxene
aug: augite
spi: spinel

"""

from sympy import solve, symbols

Mg_num_m, melt_frac = symbols('Mg_num_m melt_frac')
MgO_b, FeO_b = symbols('MgO_b FeO_b')
Xr_ol, Xr_opx, Xr_aug, Xr_spi = symbols('Xr_ol Xr_opx Xr_aug Xr_spi')
KD_olm, KD_opxm, KD_augm, KD_spim = symbols('KD_olm KD_opxm KD_augm KD_spim')

# func = solve(
#     (MgO_b - (1-melt_frac) * (
#         Xr_ol * (66.67 * Mg_num_m / (Mg_num_m + KD_olm * (1-Mg_num_m))) + 
#         Xr_opx * (50.0 * Mg_num_m / (Mg_num_m + KD_opxm * (1-Mg_num_m))) +
#         Xr_aug * (25.0) * Mg_num_m / (Mg_num_m + KD_augm * (1-Mg_num_m))) +
#         Xr_spi * (33.33 * Mg_num_m / (Mg_num_m + KD_spim * (1-Mg_num_m)))
#     ) / (
#         MgO_b + FeO_b - (1-melt_frac) * (
#             Xr_ol * 66.67 + Xr_opx * 50.0 + Xr_aug * 25.0 + Xr_spi * 33.33
#         )
#     )
# )

def _est_residual_phase_proportions(X_o: list[float], Y_o: list[float],
                                F: float) -> list[float]:
    """
    Estimates the phase proportions in the mantle residuum after primary melt extraction during
    mantle melting. This is an internal function used during iteration to calculate Mg# of the
    primary melt.
    
    This equation is from Kinzler and Grove (1992) Part 1. Appendix. Step 8:
    8. Estimate phase proportions present in solid after melt of increment F is produced (X_r) using
    the initial phase proportions (X_o) and the stoichiometric coefficients of the melting reaction
    (Y) with the following:
    X_r^ol = (X_o^ol - Y^ol * F)/(1-F)
    
    Parameters
    ----------
    X_o : list[float]
        Initial phase proportions of the original solid mantle, wt fraction.
    
    Y_o : list[float]
        Stoichiometric coefficients of the melting reaction, wt fraction.
    
    F : Melt fraction out of 1.

    Returns
    -------
    list[float]
        Estimate of the phase proportions in the mantle residuum after primary melt extraction for
        all phases, normalized, wt fraction. List order matches that as passed to X_o and Y_o.
    
    Raises
    ------
    InputError
        Raises if length of X_o does not match that of Y_o.

    Notes
    -----
    The order of the mineral values passed to X_o MUST match those passed to Y_o.
    """
    
    # check for input error
    if len(X_o) != len(Y_o):
        raise ValueError("Length of X_o must match that of Y_o.")
    
    Xr_list = []
    for xt in range(len(X_o)):
        Xr_xt = (X_o[xt] - Y_o[xt] * F) / (1-F)
        Xr_list.append(Xr_xt)
    
    sum_Xr_list = sum(Xr_list)
    Xr_list = [i/sum_Xr_list for i in Xr_list] # normalize to 1
    
    return Xr_list
    
def convert_phase_proportion_units(X_xt: list[float]) -> list[float]:
    pass

