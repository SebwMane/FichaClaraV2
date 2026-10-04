"""Pruebas de network/{weights,initialization,topology} (ANALYSIS §6.1, §6.2)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from omega.config.seeds import SeedKey, make_rng
from omega.network.initialization import random_uniform_weights
from omega.network.topology import (
    adjacency,
    binary_clustering,
    binary_degree,
    component_labels,
    component_sizes,
    giant_component_nodes,
    strength,
    submatrix,
    topology_observables,
    weighted_clustering,
)
from omega.network.weights import (
    clip_unit,
    from_upper_triangle,
    permute,
    upper_triangle,
    validate_weight_matrix,
)


def complete(n: int) -> np.ndarray:
    return np.ones((n, n)) - np.eye(n)


def two_cliques(m: int) -> np.ndarray:
    w = np.zeros((2 * m, 2 * m))
    w[:m, :m] = 1.0
    w[m:, m:] = 1.0
    np.fill_diagonal(w, 0.0)
    return w


# ---------- inicializacion ----------

def test_uniform_marginal_ks() -> None:
    w = random_uniform_weights(200, make_rng(SeedKey(2024, (0, 0, 0))))
    validate_weight_matrix(w)
    res = stats.kstest(upper_triangle(w), "uniform")
    assert res.pvalue > 0.01


def test_initialization_invariants_and_determinism() -> None:
    a = random_uniform_weights(30, make_rng(SeedKey(5, (1, 2, 3))))
    b = random_uniform_weights(30, make_rng(SeedKey(5, (1, 2, 3))))
    c = random_uniform_weights(30, make_rng(SeedKey(5, (1, 2, 4))))
    assert a.tobytes() == b.tobytes()
    assert a.tobytes() != c.tobytes()
    assert a.dtype == np.float64 and np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0.0) and a.min() >= 0.0 and a.max() < 1.0


def test_initialization_validation(rng: np.random.Generator) -> None:
    with pytest.raises(ValueError):
        random_uniform_weights(1, rng)
    with pytest.raises(TypeError):
        random_uniform_weights(5, object())  # type: ignore[arg-type]


# ---------- weights ----------

def test_validate_rejects_bad_matrices() -> None:
    good = complete(4)
    validate_weight_matrix(good)
    bad_sym = good.copy()
    bad_sym[0, 1] = 0.5
    bad_diag = good.copy()
    bad_diag[0, 0] = 0.1
    bad_range = good * 2.0
    bad_nan = good.copy()
    bad_nan[0, 1] = bad_nan[1, 0] = np.nan
    for bad in (bad_sym, bad_diag, bad_range, bad_nan, np.zeros((3, 4)), np.zeros((1, 1))):
        with pytest.raises(ValueError):
            validate_weight_matrix(bad)
    with pytest.raises(TypeError):
        validate_weight_matrix(good.astype(np.float32))


def test_upper_triangle_roundtrip_and_no_mutation(rng: np.random.Generator) -> None:
    w = random_uniform_weights(9, rng)
    before = w.copy()
    v = upper_triangle(w)
    assert v.shape == (36,)
    assert np.array_equal(from_upper_triangle(v, 9), w)
    assert np.array_equal(w, before)
    v[0] = 99.0
    assert np.array_equal(w, before)
    with pytest.raises(ValueError):
        from_upper_triangle(v, 8)


def test_clip_unit_pure() -> None:
    w = np.array([[0.0, 1.5], [-0.2, 0.0]])
    out = clip_unit(w)
    assert out.min() >= 0.0 and out.max() <= 1.0 and w[0, 1] == 1.5


def test_permute_equivariance(rng: np.random.Generator) -> None:
    w = random_uniform_weights(10, rng)
    perm = rng.permutation(10)
    pw = permute(w, perm)
    validate_weight_matrix(pw)
    assert pw[2, 5] == w[perm[2], perm[5]]
    assert np.array_equal(permute(pw, np.argsort(perm)), w)
    # observables escalares invariantes
    assert np.allclose(np.sort(strength(pw)), np.sort(strength(w)), atol=1e-12)
    assert abs(weighted_clustering(pw) - weighted_clustering(w)) <= 1e-12
    with pytest.raises(ValueError):
        permute(w, np.zeros(10, dtype=np.int64))


# ---------- topologia ----------

def test_empty_graph() -> None:
    n = 20
    o = topology_observables(np.zeros((n, n)), 0.1, 0.1)
    assert o.giant_size == 1 and o.giant_fraction == pytest.approx(1 / n)
    assert o.n_components == n and o.binary_density == 0.0
    assert o.mean_binary_degree == 0.0 and o.clustering_binary == 0.0
    assert o.clustering_weighted == 0.0 and o.cv_strength == 0.0
    assert o.frac_at_zero == 1.0 and o.frac_at_one == 0.0


def test_complete_graph() -> None:
    n = 15
    o = topology_observables(complete(n), 0.1, 0.1)
    assert o.giant_fraction == 1.0 and o.n_components == 1 and o.n_large_components == 1
    assert o.binary_density == 1.0
    assert o.clustering_binary == pytest.approx(1.0, abs=1e-12)
    assert o.clustering_weighted == pytest.approx(1.0, abs=1e-12)
    assert o.mean_binary_degree == n - 1
    assert o.mean_strength == n - 1 and o.std_strength == 0.0 and o.cv_strength == 0.0
    assert o.frac_at_one == 1.0 and o.mean_weight == 1.0


def test_two_cliques() -> None:
    w = two_cliques(10)
    o = topology_observables(w, 0.1, 0.1)
    assert o.n_large_components == 2 and o.n_components == 2
    assert o.giant_size == 10 and o.giant_fraction == 0.5
    labels = component_labels(adjacency(w, 0.1))[1]
    assert sorted(component_sizes(labels)) == [10, 10]
    giant = giant_component_nodes(labels)
    assert np.array_equal(giant, np.arange(10))  # desempate determinista: etiqueta menor
    assert submatrix(w, giant).shape == (10, 10)


def test_weighted_clustering_in_unit_interval(rng: np.random.Generator) -> None:
    for _ in range(5):
        w = random_uniform_weights(40, rng)
        c = weighted_clustering(w)
        assert 0.0 <= c <= 1.0
        # para U(0,1) i.i.d.: C_W ~ E[w]=0.5
        assert abs(c - 0.5) < 0.1


def test_binary_clustering_known_graphs() -> None:
    n = 6
    a = np.zeros((n, n), dtype=bool)
    for i in range(n):  # anillo: sin triangulos
        a[i, (i + 1) % n] = a[(i + 1) % n, i] = True
    assert binary_clustering(a) == 0.0
    tri = np.zeros((3, 3), dtype=bool)
    tri[0, 1] = tri[1, 0] = tri[1, 2] = tri[2, 1] = tri[0, 2] = tri[2, 0] = True
    assert binary_clustering(tri) == 1.0
    assert list(binary_degree(tri)) == [2, 2, 2]


def test_adjacency_threshold_and_no_mutation() -> None:
    w = np.array([[0.0, 0.1, 0.05], [0.1, 0.0, 0.5], [0.05, 0.5, 0.0]])
    before = w.copy()
    a = adjacency(w, 0.1)
    assert a[0, 1] and a[1, 2] and not a[0, 2] and not a.diagonal().any()
    assert np.array_equal(w, before)
    with pytest.raises(ValueError):
        adjacency(w, 1.5)


def test_topology_observables_random_in_ranges(rng: np.random.Generator) -> None:
    w = random_uniform_weights(100, rng)
    o = topology_observables(w, 0.1, 0.1)
    assert o.n == 100 and o.giant_fraction == 1.0
    assert 0.0 <= o.binary_density <= 1.0
    assert o.binary_density == pytest.approx(0.9, abs=0.02)
    assert o.cv_strength < 0.2
    assert o.frac_at_zero == 0.0
