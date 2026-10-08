import numpy as np
from scipy import sparse

from omega.c0.references import rng_from_key
from omega.diagnostics.ball_growth import ball_profile
from omega.diagnostics.sampled_growth import (
    dimension_status,
    giant_component,
    report_category,
    sampled_ball_profile,
    spectral_return,
)


def _ring(n: int, k: int) -> sparse.csr_array:
    i = np.arange(n)
    r = np.concatenate([i for _ in range(k // 2)])
    c = np.concatenate([(i + d) % n for d in range(1, k // 2 + 1)])
    return sparse.csr_array(sparse.coo_array((np.ones(r.size), (r, c)), shape=(n, n)))


def _torus(side: int) -> sparse.csr_array:
    n = side * side
    v = np.arange(n)
    i, j = v // side, v % side
    r = np.concatenate([v, v])
    c = np.concatenate([((i + 1) % side) * side + j, i * side + (j + 1) % side])
    return sparse.csr_array(sparse.coo_array((np.ones(r.size), (r, c)), shape=(n, n)))


def test_ring_converges_class_1() -> None:
    p = sampled_ball_profile(_ring(2000, 12), rng_from_key((1, 1)), n_sources=50)
    assert p["level1"]["status"] == "CV_CERO"
    assert p["level2"]["0.05"]["status"] == "CONVERGE"
    assert p["level2"]["0.05"]["class"] == 1
    assert report_category(p["level1"]["pass"], "CONVERGE", 1) == "GEOMETRIA_GRUESA(1)"


def test_square_torus_class_2() -> None:
    p = sampled_ball_profile(_torus(60), rng_from_key((1, 2)), n_sources=60)
    assert p["level1"]["status"] in ("CV_CERO", "HOMOGENEIZA")
    assert p["level2"]["0.05"]["class"] == 2


def test_random_regular_has_no_window() -> None:
    import networkx as nx

    g = nx.random_regular_graph(12, 2000, seed=3)
    a = sparse.csr_array(nx.to_scipy_sparse_array(g))
    p = sampled_ball_profile(a, rng_from_key((1, 3)), n_sources=50)
    assert p["level1"]["status"] == "SIN_VENTANA"
    assert report_category(False, p["level2"]["0.05"]["status"], None) == "NO_GEOMETRICO"


def test_rising_dimension_is_cruce() -> None:
    d = [1.0 + 0.2 * i for i in range(8)]
    assert dimension_status(d, 0.05)["status"] == "CRUCE"
    assert dimension_status(d[:3], 0.05)["status"] == "SIN_VENTANA_D"
    flat = [2.0, 2.3, 2.1, 2.05, 2.03, 2.02]
    assert dimension_status(flat, 0.05)["status"] == "CONVERGE"


def test_sampled_matches_exact_when_all_sources() -> None:
    from omega.controls.random_geometric import rgg_torus

    dense = rgg_torus(300, 2, 10.0, rng_from_key((5, 5)))
    exact = ball_profile(dense)
    g = giant_component(sparse.csr_array(dense))
    assert g.shape[0] == exact["n_giant"]
    p = sampled_ball_profile(sparse.csr_array(dense), rng_from_key((5, 6)), n_sources=g.shape[0])
    k = min(len(exact["mean"]), len(p["mean"]), 6)
    assert np.allclose(p["mean"][:k], exact["mean"][:k])
    assert np.allclose(p["cv"][:k], exact["cv"][:k])


def test_spectral_ring_dimension_near_1() -> None:
    s = spectral_return(_ring(600, 4), rng_from_key((2, 2)), n_sources=8, t_max=300)
    assert s["D_s_median_upper"] is not None
    assert 0.6 < s["D_s_median_upper"] < 1.4
