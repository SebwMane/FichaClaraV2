"""Tests unitarios del clasificador estructural R2.4 (omega.landscape.structure)."""

from __future__ import annotations

import networkx as nx
import numpy as np

from omega.experiments.reference_graphs import periodic_lattice
from omega.landscape.kkt import codegree, open_p3_count
from omega.landscape.structure import (
    CLASS_EMPTY,
    CLASS_MULTI_CLIQUE,
    CLASS_OTHER,
    CLASS_OVERLAPPING,
    CLASS_SINGLE_CLIQUE,
    CLASS_UNIFORM,
    classify_structure,
    kappa_inj,
)
from omega.types import FloatArray


def _blocks(n: int, groups: list[range]) -> FloatArray:
    w = np.zeros((n, n))
    for g in groups:
        idx = np.array(list(g))
        w[np.ix_(idx, idx)] = 1.0
    np.fill_diagonal(w, 0.0)
    return w


def test_empty() -> None:
    assert classify_structure(np.zeros((10, 10))).state_class == CLASS_EMPTY
    w = np.full((10, 10), 5e-7)
    np.fill_diagonal(w, 0.0)
    assert classify_structure(w).state_class == CLASS_EMPTY


def test_uniform() -> None:
    w = np.full((12, 12), 0.3)
    np.fill_diagonal(w, 0.0)
    c = classify_structure(w)
    assert c.state_class == CLASS_UNIFORM
    assert not c.has_halo or c.halo_is_clique_union  # uniforme 0.3: estrato debil completo


def test_single_clique() -> None:
    c = classify_structure(_blocks(30, [range(0, 10)]))
    assert c.state_class == CLASS_SINGLE_CLIQUE
    assert c.component_sizes == (10,)
    assert not c.has_halo and c.halo_is_clique_union is None
    assert abs(c.kappa_inj[0] - 1.0) < 1e-12
    assert c.open_p3 == 0


def test_multi_clique() -> None:
    c = classify_structure(_blocks(40, [range(0, 10), range(10, 18), range(18, 30)]))
    assert c.state_class == CLASS_MULTI_CLIQUE
    assert c.component_sizes == (10, 8, 12)
    assert all(abs(k - 1.0) < 1e-12 for k in c.kappa_inj)
    assert c.open_p3 == 0


def test_two_cliques_overlapping_in_two_vertices() -> None:
    w = _blocks(50, [range(0, 20), range(18, 44)])  # tamanos 20 y 26, comparten los vertices 18 y 19
    c = classify_structure(w)
    assert c.state_class == CLASS_OVERLAPPING
    assert c.component_sizes == (44,)
    assert c.components[0].kind == "solapada" and c.components[0].n_big_cliques == 2
    assert c.kappa_inj[0] < 1.0
    assert c.open_p3 > 0  # el soporte no es union disjunta de cliques


def test_clique_plus_halo() -> None:
    n = 100
    w = np.zeros((n, n))
    w[22:, 22:] = 0.0879
    w[:22, :22] = 1.0
    np.fill_diagonal(w, 0.0)
    c = classify_structure(w)
    assert c.state_class == CLASS_SINGLE_CLIQUE
    assert c.has_halo and c.halo_is_clique_union is True
    assert c.component_sizes == (22,)
    assert c.open_p3 == 0


def test_halo_not_clique_union() -> None:
    w = _blocks(20, [range(0, 6)])
    w[10, 11] = w[11, 10] = 0.3
    w[11, 12] = w[12, 11] = 0.3  # P3 abierto 10-11-12 en el estrato debil
    c = classify_structure(w)
    assert c.state_class == CLASS_SINGLE_CLIQUE
    assert c.has_halo and c.halo_is_clique_union is False
    assert c.open_p3 == 1


def test_torus_6cubed_is_other() -> None:
    w = periodic_lattice((6, 6, 6))
    c = classify_structure(w)
    assert c.state_class == CLASS_OTHER
    assert c.component_sizes == (216,)
    assert c.components[0].kind == "otro"
    assert c.open_p3 > 0


def test_erdos_renyi_is_other() -> None:
    g = nx.gnp_random_graph(60, 0.1, seed=7)
    w = np.asarray(nx.to_numpy_array(g), dtype=np.float64)
    assert classify_structure(w).state_class == CLASS_OTHER


def test_weak_only_non_uniform_is_other() -> None:
    rng = np.random.Generator(np.random.PCG64(1))
    a = rng.random((15, 15)) * 0.4 + 0.01
    w = np.triu(a, 1)
    w = w + w.T
    assert classify_structure(w).state_class == CLASS_OTHER


def test_kappa_inj_and_helpers() -> None:
    assert abs(kappa_inj(np.ones((7, 7), dtype=bool) & ~np.eye(7, dtype=bool)) - 1.0) < 1e-12
    cyc = np.asarray(nx.to_numpy_array(nx.cycle_graph(6)), dtype=bool)
    assert kappa_inj(cyc) == 0.0
    w = periodic_lattice((3, 3, 3))
    assert codegree(w).shape == w.shape and np.all(np.diag(codegree(w)) == 0.0)
    assert open_p3_count(np.ones((5, 5)) - np.eye(5), 0.5) == 0
