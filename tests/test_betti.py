"""Tests de omega.topology.betti (WP-C)."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pytest

from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.topology.betti import (
    clique_complex_betti,
    count_four_cycles,
    count_triangles,
    gf2_rank,
    short_cycle_betti1,
)


def _adj(n: int, edges: list[tuple[int, int]]) -> np.ndarray:
    a = np.zeros((n, n), dtype=np.bool_)
    for i, j in edges:
        a[i, j] = a[j, i] = True
    return a


def test_gf2_rank() -> None:
    assert gf2_rank([]) == 0
    assert gf2_rank([0b011, 0b110, 0b101]) == 2
    assert gf2_rank([1 << 200, 1 << 100, (1 << 200) | (1 << 100)]) == 2  # enteros grandes


def test_clique_complex_k4_c5_octahedron_disjoint() -> None:
    r = clique_complex_betti(complete_graph(4) > 0)
    assert r.betti == (1, 0, 0) and r.status == "ok" and r.counts == (4, 6, 4, 1)
    r = clique_complex_betti(_adj(5, [(i, (i + 1) % 5) for i in range(5)]))
    assert r.betti == (1, 1, 0)
    octa = np.ones((6, 6), dtype=np.bool_)
    np.fill_diagonal(octa, False)
    for i, j in [(0, 1), (2, 3), (4, 5)]:
        octa[i, j] = octa[j, i] = False
    r = clique_complex_betti(octa)
    assert r.betti == (1, 0, 1) and r.euler_betti == 2 and r.euler_counts == 2
    two = _adj(6, [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    r = clique_complex_betti(two)
    assert r.betti[0] == 2 and r.betti[1] == 0


def test_clique_over_budget_and_validation() -> None:
    r = clique_complex_betti(complete_graph(60) > 0, 2, 100)
    assert r.status == "over_budget"
    with pytest.raises(ValueError):
        clique_complex_betti(np.array([[False, True], [False, False]]))


def test_counts_closed_forms_match_enumeration() -> None:
    rng = np.random.Generator(np.random.PCG64(3))
    a = np.triu(rng.random((14, 14)) < 0.35, 1)
    a = a | a.T
    n = a.shape[0]
    tri = sum(1 for i, j, k in combinations(range(n), 3) if a[i, j] and a[j, k] and a[i, k])
    assert count_triangles(a) == tri
    c4 = 0
    for q in combinations(range(n), 4):
        i, j, k, l = q
        for cyc in ((i, j, k, l), (i, j, l, k), (i, k, j, l)):
            if all(a[cyc[t], cyc[(t + 1) % 4]] for t in range(4)):
                c4 += 1
    assert count_four_cycles(a) == c4
    assert count_four_cycles(complete_graph(4) > 0) == 3


def test_short_cycle_ring_torus2_complete() -> None:
    assert short_cycle_betti1(periodic_lattice((9,)) > 0).b1 == 1
    r = short_cycle_betti1(periodic_lattice((6, 6)) > 0)
    assert r.b1 == 2 and r.b0 == 1 and r.n_edges == 72 and r.status == "ok"
    assert r.b1_density == pytest.approx(2 / 72)
    for n in (4, 5, 8):
        assert short_cycle_betti1(complete_graph(n) > 0).b1 == 0


def test_short_cycle_torus3() -> None:
    assert short_cycle_betti1(periodic_lattice((6, 6, 6)) > 0).b1 == 3


def test_short_cycle_small_torus_measured_deviation() -> None:
    # Medido: en 4^3 (y 4^2) el lazo de envoltura tiene longitud 4 y queda relleno: beta1 = 0, no 3.
    # El diseno pide 3; se informa como desviacion (tamaño finito) sin ensanchar.
    assert short_cycle_betti1(periodic_lattice((4, 4, 4)) > 0).b1 == 0
    assert short_cycle_betti1(periodic_lattice((5, 5, 5)) > 0).b1 == 3


@pytest.mark.slow
def test_short_cycle_torus3_9() -> None:
    assert short_cycle_betti1(periodic_lattice((9, 9, 9)) > 0).b1 == 3


def test_short_cycle_triangles_only_and_budget() -> None:
    sq = periodic_lattice((4,)) > 0  # C4
    assert short_cycle_betti1(sq, 3).b1 == 1
    assert short_cycle_betti1(sq, 4).b1 == 0
    r = short_cycle_betti1(complete_graph(30) > 0, 4, 100)
    assert r.status == "over_budget" and r.n_faces > 100
    with pytest.raises(ValueError):
        short_cycle_betti1(sq, 5)


def test_short_cycle_empty_graph() -> None:
    r = short_cycle_betti1(np.zeros((5, 5), dtype=np.bool_))
    assert (r.b0, r.b1, r.n_edges, r.b1_density) == (5, 0, 0, 0.0)
