import pytest
from stellar_geology import kinzler_grove

def test_estimate_of_phase_proportions_after_melting_step():
    """ Tests calculation against that given in the example in the Appendix of Kinzler and Grove
    (1992) part 1.
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