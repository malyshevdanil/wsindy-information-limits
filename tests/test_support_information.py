import numpy as np
from src.support_information import efficient_column_information


def test_projection_information_matches_schur_complement():
    rng = np.random.default_rng(313); W = rng.normal(size=(100, 7))
    I = W.T@W; nuisance = [0, 1, 2, 4, 5, 6]
    exact = I[3, 3]-I[3, nuisance]@np.linalg.solve(I[np.ix_(nuisance, nuisance)], I[nuisance, 3])
    result = efficient_column_information(W, 3, nuisance, 1e-10)
    np.testing.assert_allclose(result['information'], exact, rtol=1e-13)


def test_full_dictionary_sensitivity_contains_original_vdp_columns():
    from src.support_information import augmented_sensitivities
    from src.path_scoring import vdp_sensitivities
    t = np.linspace(0, .2, 9); initial = [2., 0.]
    x, S = augmented_sensitivities(t, initial)
    original, U = vdp_sensitivities(t, [1., -1., 1., -1.], initial)
    np.testing.assert_allclose(x, original, atol=1e-9)
    np.testing.assert_allclose(S[:, :, [1, 9, 10, 15, 18, 19]], U, atol=1e-9)
