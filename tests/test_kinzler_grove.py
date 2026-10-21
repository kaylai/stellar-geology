import pytest
from stellar_geology import kinzler_grove

def test_estimate_of_phase_proportions_after_melting_step():
    """ Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 1, step 8.
    """
    # published values
    initial_phase_prop = [
        0.193, # augite
        0.538, # olivine
        0.239, # opx
        0.029, # spinel
    ]
    
    stoich_coefficients = [
        0.97, # augite
        -0.30, # olivine
        0.25, # opx
        0.08, # spinel
    ]
    
    melt_frac = 0.1
    
    expected_res_phase_prop = [
        0.107, # augite
        0.631, # olivine
        0.239, # opx
        0.023, # spinel
    ]
    
    res_phase_prop = kinzler_grove._est_residual_phase_proportions(
        X_o=initial_phase_prop,
        Y_o=stoich_coefficients,
        F=melt_frac)
    
    for i in range(len(res_phase_prop)):
        assert (
            res_phase_prop[i] == pytest.approx(expected_res_phase_prop[i], rel=1e-1)
        ), f"\n res_phase_prop[{i}] was {res_phase_prop[i]}, expected {expected_res_phase_prop[i]}"

def test_estimate_liquid_TiO2_and_alkalis():
    """Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 2, step 1. 
    """
    # published values
    C_L_pub = {
        "TiO2": 0.854,
        "Na2O": 2.185,
        "K2O": 0.079,
    }
    
    # inputs from Kinzler and Grove (1992) part 2
    C_o_vals = {
        "TiO2": 0.174,
        "Na2O": 0.281,
        "K2O": 0.008,
    }
    
    D_B_vals = {
            "TiO2": 0.147,
            "Na2O": 0.049,
            "K2O": 0.001,
        }
    
    P_B_vals = {
            "TiO2": 0.433,
            "Na2O": 0.204,
            "K2O": 0.001,
        }
    
    F_val = 0.1
    
    # calculate new values
    C_L_calc = kinzler_grove._est_liquid_TiO2_and_alkalis(C_o=C_o_vals, D_B=D_B_vals, P_B=P_B_vals,
                                                          F=F_val)
    
    for k in C_L_calc.keys():
        assert (
            C_L_calc[k] == pytest.approx(C_L_pub[k], rel=1e-2)
        ), f"\n The liquid concentration of {k} was {C_L_calc[k]}, expected {C_L_pub[k]}"

@pytest.mark.parametrize("missing", ["TiO2", "Na2O", "K2O"])
def test_missing_keys_raises(missing):
    required_keys = {"TiO2": 1.0, "Na2O": 2.0, "K2O": 0.5}
    missing_keys = {k: v for k, v in required_keys.items() if k != missing}
    with pytest.raises(ValueError, match=missing):
        kinzler_grove._est_liquid_TiO2_and_alkalis(C_o=missing_keys, D_B=missing_keys,
                                                   P_B=missing_keys, F=1)

def test_extra_keys_warns():
    required_keys = {"TiO2": 1.0, "Na2O": 2.0, "K2O": 0.5}
    extra = {**required_keys, "SiO2": 50.0}
    with pytest.warns(UserWarning, match="SiO2"):
        kinzler_grove._est_liquid_TiO2_and_alkalis(C_o=extra, D_B=extra, P_B=extra, F=1)

def test_required_keys_only_no_warning():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes a failure
        required_keys = {"TiO2": 1.0, "Na2O": 2.0, "K2O": 0.5}
        kinzler_grove._est_liquid_TiO2_and_alkalis(C_o=required_keys, D_B=required_keys,
                                                   P_B=required_keys, F=1)

def test_estimate_D_CaO_aug_opx():
    """Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 2, step 1. 
    """
    # published values
    D_CaO_aug_pub = 1.267
    D_CaO_opx_pub = 0.204
    
    # inputs from Kinzler and Grove (1992) part 2
    C_L_Na2O = 2.185
    
    D_CaO_aug_calc, D_CaO_opx_calc = kinzler_grove._est_D_CaO_aug_opx(C_L_Na2O)
    
    assert D_CaO_aug_calc == pytest.approx(D_CaO_aug_pub, rel=1e-2)
    assert D_CaO_opx_calc == pytest.approx(D_CaO_opx_pub, rel=1e-2)
    
def test_estimate_D_B_and_P_B_CaO():
    """Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 2, step 1. 
    """
    # published values
    D_B_CaO_pub = 0.310
    P_B_CaO_pub = 1.272
    
    # inputs from Kinzler and Grove (1992) part 2
    D_CaO_aug_val = 1.267
    D_CaO_opx_val = 0.204
    X_vals = {
        "olivine": 0.538,
        "orthopyroxene": 0.24,
        "augite": 0.193,
        "spinel": 0.029
    }
    Y_vals = {
        "olivine": -0.3,
        "orthopyroxene": 0.25,
        "augite": 0.97,
        "spinel": 0.08
    }
    
    D_B_CaO, P_B_CaO = kinzler_grove._est_D_B_and_P_B_CaO(D_CaO_aug_val, D_CaO_opx_val,
                                                          X_vals, Y_vals)
    
    assert D_B_CaO == pytest.approx(D_B_CaO_pub, rel=1e-2)
    assert P_B_CaO == pytest.approx(P_B_CaO_pub, rel=1e-2)
    
def test_estimate_C_L_CaO():
    """Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 2, step 1.
    """
    # published values
    C_L_CaO_pub = 11.287
    
    # inputs from Kinzler and Grove (1992) part 2
    C_o_CaO_val = 3.192
    D_B_CaO_val = 0.310
    P_B_CaO_val = 1.272
    F_val = 0.1
    
    C_L_CaO_calc = kinzler_grove._est_C_L_CaO(C_o_CaO_val, D_B_CaO_val, P_B_CaO_val, F_val)
    
    assert C_L_CaO_calc == pytest.approx(C_L_CaO_pub, rel=1e-2)

def test_estimate_NaKnum():
    """Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 2, step 1.
    """
    # published values
    NaKnum_L = 0.167
    
    # inputs from Kinzler and Grove (1992) part 2
    C_L_Na2O_val = 2.185
    C_L_K2O_val = 0.079
    C_L_CaO_val = 11.287
    
    NaKnum_calc = kinzler_grove._est_NaKnum(C_L_Na2O_val, C_L_K2O_val, C_L_CaO_val)