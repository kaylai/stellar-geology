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
import warnings

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

def _est_liquid_TiO2_and_alkalis(C_o: dict[str, float], D_B: dict[str, float],
                                 P_B: dict[str, float], F:float) -> dict[str, float]:
    """
    Estimates the wt% concentration of TiO2, Na2O, and K2O in the liquid of a primary mantle melt
    using the nonmodal batch melting equation.
    
    This equation is from Kinzler and Grove (1992) Part 2. Appendix. Step 1:
    
    C_i,L = C_o / (D_i,B + F*(1-P_i,B))
    
    where C_i,L is the concentration of oxide i in the melt; C_o is the initial concentration of
    the oxide in the bulk solid; D_B is the bulk partition coefficient for the oxide, weighted by
    the proportion of minerals in the initial solid; F is the fraction of liquid; P_B is the bulk
    partition coefficient of the oxide, weighted by proportion of minerals going into the melt phase.
    
    Parameters
    ----------
    C_o: dict[str, float]
        Initial concentrations of oxides TiO2, Na2O, and K2O, wt%.
    
    D_B: dict[str, float]
        Bulk partition coefficients for oxides TiO2, Na2O, and K2O weighted by the proportion of
        minerals in the initial solid.
        
    P_B: dict[str, float]
        Bulk partition coefficients for oxides TiO2, Na2O, and K2O weighted by the proportion of
        minerals going into the melt phase.
    
    F: float
        Fraction of liquid produced, out of 1.
    
    Returns
    -------
    dict[str, float]
        Estimate of the wt% concentrations of TiO2, Na2O, and K2O in a primary mantle melt.
    
    Raises
    ------
    ValueError
        Raised if any of key "TiO2", "Na2O", or "K2O" is missing from any argument.
    
    Warns
    -----
    UserWarning
        Warned if any key is in any argument in addition to the required keys "TiO2", "Na2O", and
        "K2O".
    
    Notes
    -----
    """
    required_keys = {"TiO2", "Na2O", "K2O"}
    def _check_keys(arg, name):
        missing = required_keys - arg.keys()
        if missing:
            raise ValueError(f"{name} is missing {sorted(missing)}")
        
        extra = arg.keys() - required_keys
        if extra:
            warnings.warn(f"{name} ignoring {sorted(extra)}", stacklevel=3)
    
    for name, arg in {"C_o": C_o, "D_B": D_B, "P_B": P_B}.items():
        _check_keys(arg, name)
    
    C_L = {}
    for ox in required_keys:
        C_L[ox] = C_o[ox] / (D_B[ox]+F*(1-P_B[ox]))
    
    print(C_L)
    return C_L

def _est_D_CaO_aug_opx(C_L_Na2O: float) -> tuple[float, float]:
    """
    Estimates the bulk partition coefficient D_CaO for augite and opx as a function of Na2O
    concentration in the melt.
    
    These equations are from Kinzler and Grove (1992) Part 2. Appendix. Step 1:
    
    D_CaO_aug = 1.07 + 0.09 * C_L_Na2O
    D_CaO_opx = 0.16 + 0.02 * C_L_Na2O
    
    Parameters
    ----------
    C_L_Na2O: float
        Concentration of Na2O in the primary liquid, in wt%
    
    Returns
    -------
    D_CaO_aug: float
        Partition coefficient of CaO for augite
    D_CaO_opx: float
        Partition coefficient of CaO for orthopyroxene
    """
    D_CaO_aug = 1.07 + 0.09 * C_L_Na2O
    D_CaO_opx = 0.16 + 0.02 * C_L_Na2O
    
    return D_CaO_aug, D_CaO_opx

def _est_D_B_and_P_B_CaO(D_CaO_aug: float, D_CaO_opx: float,
                         X_i: dict[str,float], Y_i: dict[str, float]) -> tuple[float, float]:
    """
    Estimates the weighted partition coefficiencts, D_B and P_B for CaO.
    
    This equation is from Kinzler and Grove (1992) Part 2. Appendix. Step 1:
    
    D_B_CaO = D_CaO_aug*X_aug + D_CaO_opx*X_opx + X_ol*0.03 + X_spi*0.029
    P_B_CaO = D_CaO_aug*Y_aug + D_CaO_opx*Y_opx + Y_ol*0.03 + Y_spi*0.018
    
    Parameters
    ----------
    D_CaO_aug: float
        Partition coefficient of CaO for augite
        
    D_CaO_opx: float
        Partition coefficient of CaO for orthopyroxene
        
    X_i: dict[str, float]
        Phase proportions of olivine, orthopyroxene, augite, and spinel in wt. fraction. Keys must
        be: "olivine", "orthopyroxene", "olivine", and "spinel".
    
    Y_i: dict[str, float]
        Melt reaction coefficients for olivine, orthopyroxene, autire, and spinel. Keys must be
        "olivine", "orthopyroxene", "olivine", and "spinel".
    
    Returns
    -------
    D_B_CaO: float
        Bulk partition coefficient fo CaO weighted by the proportion of minerals in the initial
        solid
    P_B_CaO: float
        Bulk partition coefficient of CaO weighted by the proportion of minerals going into the
        melt phase
    """
    D_B_CaO = (D_CaO_aug*X_i["augite"] + D_CaO_opx*X_i["orthopyroxene"] + 
               X_i["olivine"]*0.03 + X_i["spinel"]*0.029)
    
    P_B_CaO = (D_CaO_aug*Y_i["augite"] + D_CaO_opx*Y_i["orthopyroxene"] +
               Y_i["olivine"]*0.03 + Y_i["spinel"]*0.018)
    
    return (D_B_CaO, P_B_CaO)

def _est_C_L_CaO(C_o_CaO: float, D_B_CaO: float, P_B_CaO: float, F: float) -> float:
    """Estimates the concentration of CaO in the primary melt.
    This equation is from Kinzler and Grove (1992) Part 2. Appendix. Step 1:
    
    C_L_CaO = C_o_CaO/(D_B_CaO+F*(1-P_B_CaO))
    
    Parameters
    ----------
    C_o_CaO: float
        Concentration of CaO in the original bulk solid
    
    D_B_CaO: float
        Bulk partition coefficient fo CaO weighted by the proportion of minerals in the initial
        solid
        
    P_B_CaO: float
        Bulk partition coefficient of CaO weighted by the proportion of minerals going into the
        melt phase
    
    F: float
        Fraction of liquid produced, out of 1.
    
    Returns
    -------
    float
        Concentration of CaO in the primary melt.
    """
    return C_o_CaO / (D_B_CaO + F*(1-P_B_CaO))

def _est_NaKnum(C_L_Na2O: float, C_L_K2O: float, C_L_CaO: float):
    """Estimates the NaK# in the primary melt.
    
    This equation is from Kinzler and Grove (1992) Part 2. Appendix. Step 1:
    
    NaK# = (Na2O+K2O)/(Na2O+K2O+CaO)
    
    Parameters
    ----------
    C_L_Na2O: float
        Concentration of Na2O in the primary melt.
    
    C_L_K2O: float
        Concentration of K2O in the primary melt.
    
    C_L_CaO: float
        Concentration of CaO in the primary melt.
    """
    return (C_L_Na2O+C_L_K2O) / (C_L_Na2O+C_L_K2O+C_L_CaO)
    
    

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
    ValueError
        Raised if length of X_o does not match that of Y_o.

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

