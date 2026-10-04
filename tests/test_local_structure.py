"""Pruebas de omega.geometry.local_structure (DESIGN §1.4-§1.5, §3.2 WP-B)."""

from __future__ import annotations

import math

import networkx as nx
import numpy as np
import pytest

from omega.config.settings import DimensionConfig
from omega.config.settings11 import LocalStructureConfig
from omega.experiments.reference_graphs import (
    complete_graph,
    periodic_lattice,
    random_geometric_torus,
)
from omega.geometry.dimension import effective_dimension
from omega.geometry.distances import hop_distance_matrix
from omega.geometry.local_structure import (
    annulus_connectivity,
    ball_mds_ratio,
    edge_detour_fraction,
    homogeneity,
    isotropy,
    locality,
    node_dimensions,
)
from omega.types import BoolArray, DimensionEstimate, FloatArray

CFG = LocalStructureConfig()
DIM = DimensionConfig()


def _setup(w: FloatArray) -> tuple[BoolArray, FloatArray, DimensionEstimate]:
    a = np.asarray(w > 0.5, dtype=np.bool_)
    d = hop_distance_matrix(a)
    return a, d, effective_dimension(d, DIM)


def _tri_lattice(n: int) -> BoolArray:
    g = nx.triangular_lattice_graph(n, n, periodic=False)
    return np.asarray(nx.to_numpy_array(g, dtype=np.float64) > 0, dtype=np.bool_)


def _tri_torus(n: int) -> BoolArray:
    a = np.zeros((n * n, n * n), dtype=np.bool_)
    for i in range(n):
        for j in range(n):
            u = i * n + j
            for di, dj in ((1, 0), (0, 1), (1, 1)):
                v = ((i + di) % n) * n + (j + dj) % n
                a[u, v] = a[v, u] = True
    return a


@pytest.fixture(scope="module")
def cube9() -> tuple[BoolArray, FloatArray, DimensionEstimate]:
    return _setup(periodic_lattice((9, 9, 9)))


def test_detour_ring_zero_and_dense_one(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    for n in (5, 6, 30):
        assert edge_detour_fraction(periodic_lattice((n,)) > 0) == 0.0 if n >= 5 else True
    assert edge_detour_fraction(cube9[0]) == 1.0
    assert edge_detour_fraction(complete_graph(8) > 0) == 1.0
    assert edge_detour_fraction(_tri_torus(6)) == 1.0
    assert edge_detour_fraction(_tri_lattice(5)) == 1.0
    assert math.isnan(edge_detour_fraction(np.zeros((4, 4), dtype=np.bool_)))


def test_detour_formula_matches_bfs() -> None:
    for seed in range(5):
        g = nx.gnp_random_graph(14, 0.25, seed=seed)
        a = np.asarray(nx.to_numpy_array(g, nodelist=range(14)) > 0, dtype=np.bool_)
        edges = list(g.edges())
        if not edges:
            continue
        good = 0
        for u, v in edges:
            h = g.copy()
            h.remove_edge(u, v)
            try:
                good += int(nx.shortest_path_length(h, u, v) <= 3)
            except nx.NetworkXNoPath:
                pass
        assert edge_detour_fraction(a) == pytest.approx(good / len(edges))


def test_locality_report() -> None:
    assert locality(periodic_lattice((20,)) > 0, CFG).ok is False
    rep = locality(periodic_lattice((6, 6)) > 0, CFG)
    assert rep.ok and rep.detour_fraction == 1.0
    assert not locality(np.zeros((3, 3), dtype=np.bool_), CFG).ok


def test_homogeneity_torus_zero_cv(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    _, d, est = cube9
    rep = homogeneity(d, est, CFG)
    assert rep.ok and rep.cv <= 1e-12 and rep.n_valid == 729
    assert rep.mu == pytest.approx(est.value, abs=0.2)


def test_homogeneity_nan_and_unusable_estimate(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    _, d, est = cube9
    di = node_dimensions(d, np.array([1.0, 2.0, 3.0, 20.0]))  # S(20)=0 => NaN en todos
    assert np.all(np.isnan(di))
    bad = effective_dimension(hop_distance_matrix(complete_graph(30) > 0), DIM)
    rep = homogeneity(d, bad, CFG)
    assert not rep.ok and rep.n_valid == 0


def test_isotropy_torus_and_slab(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    _, d, est = cube9
    rep = isotropy(d, est, 3, CFG, np.random.Generator(np.random.PCG64(0)))
    assert rep.ok and rep.median_ratio >= 0.99 and rep.n_sources == 64 and rep.radius == 5
    _, ds, es = _setup(periodic_lattice((32, 5, 5)))
    slab = isotropy(ds, es, 3, CFG, np.random.Generator(np.random.PCG64(0)))
    assert not slab.ok
    assert slab.median_ratio == pytest.approx(0.27, abs=0.02)


def test_isotropy_class_one_and_none(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    _, d, est = cube9
    assert isotropy(d, est, 1, CFG, np.random.Generator(np.random.PCG64(0))).ok
    assert not isotropy(d, est, None, CFG, np.random.Generator(np.random.PCG64(0))).ok


def test_isotropy_deterministic_by_rng() -> None:
    _, d, est = _setup(periodic_lattice((9, 9, 9)))
    r1 = isotropy(d, est, 3, CFG, np.random.Generator(np.random.PCG64(7)))
    r2 = isotropy(d, est, 3, CFG, np.random.Generator(np.random.PCG64(7)))
    assert r1 == r2


def test_ball_mds_ratio_cases(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    _, d, _ = cube9
    assert ball_mds_ratio(d, 0, 4, 1) == pytest.approx(1.0)
    assert ball_mds_ratio(d, 0, 4, 3) == pytest.approx(1.0, abs=1e-9)
    assert ball_mds_ratio(d, 0, 0, 3) == 0.0  # bola de un solo nodo
    ring = hop_distance_matrix(periodic_lattice((200,)) > 0)
    assert ball_mds_ratio(ring, 0, 10, 2) < 0.2   # bola de un anillo: casi unidimensional
    with pytest.raises(ValueError):
        ball_mds_ratio(d, -1, 3, 2)


def test_annulus_torus_ok_and_tree_fails(cube9: tuple[BoolArray, FloatArray, DimensionEstimate]) -> None:
    a, d, est = cube9
    rep = annulus_connectivity(a, d, est, CFG, np.random.Generator(np.random.PCG64(0)))
    assert rep.ok and min(rep.fractions.values()) == 1.0 and set(rep.fractions) == {2, 3, 4}
    g = nx.balanced_tree(3, 5)
    at, dt, et = _setup(nx.to_numpy_array(g))
    rt = annulus_connectivity(at, dt, et, CFG, np.random.Generator(np.random.PCG64(0)))
    assert not rt.ok and max(rt.fractions.values()) == 0.0


def test_annulus_scale_insufficient_and_ring() -> None:
    a, d, est = _setup(periodic_lattice((300,)))
    rep = annulus_connectivity(a, d, est, CFG, np.random.Generator(np.random.PCG64(0)))
    assert not rep.ok and all(v == 0.0 for v in rep.fractions.values())
    small = DimensionEstimate(2.0, 0.0, "ok", (0, 0), np.array([1.0]), np.zeros(1), np.zeros(1), True, "shell")
    rep2 = annulus_connectivity(a, d, small, CFG, np.random.Generator(np.random.PCG64(0)))
    assert not rep2.ok and rep2.fractions == {}


@pytest.mark.slow
def test_rgg3_reference_values() -> None:
    w = random_geometric_torus(800, 3, 12, np.random.Generator(np.random.PCG64(0)), euclidean=False)
    a, d, est = _setup(w)
    assert homogeneity(d, est, CFG).cv < 0.2
    assert locality(a, CFG).ok
    assert isotropy(d, est, 3, CFG, np.random.Generator(np.random.PCG64(0))).ok
    assert annulus_connectivity(a, d, est, CFG, np.random.Generator(np.random.PCG64(0))).ok
