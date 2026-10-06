"""Omega-D L-OD-0: solapamiento de difusion (cotas, simetria, degeneracion de cliques, referencias)."""

from __future__ import annotations

from typing import Any

import numpy as np

from omega.c0 import references as R
from omega.diffusion.overlap import edge_scaling, giant_dense, lazy_sym_operator, overlap_powers


def _measure(adj: np.ndarray) -> tuple[dict[int, np.ndarray], dict[str, Any]]:
    a = giant_dense(adj)
    al = overlap_powers(lazy_sym_operator(a))
    iu, ju = np.nonzero(np.triu(a, 1))
    return al, edge_scaling(al, iu, ju)["summary"]


def _square(side: int) -> np.ndarray:
    n = side * side
    a = np.zeros((n, n))
    for i in range(side):
        for j in range(side):
            for u in (((i + 1) % side) * side + j, i * side + (j + 1) % side):
                a[i * side + j, u] = a[u, i * side + j] = 1.0
    return a


def test_bounds_and_unit_diagonal() -> None:
    al, _ = _measure(R.random_regular(200, 6, R.rng_from_key((1, 2, 3))))
    for m in al.values():
        assert m.min() >= 0.0 and m.max() <= 1.0
        assert np.allclose(m, m.T)
        assert np.allclose(np.diag(m), 1.0)


def test_lattices_scaling_values_reported() -> None:
    _, s1 = _measure(R.ring_lattice(400, 2))
    _, s2 = _measure(_square(30))
    print("\nring  Q_8 median", s1["median_Q"]["8"], "s", s1["median_s"], "(d/2 = 0.5)")
    print("square Q_8 median", s2["median_Q"]["8"], "s", s2["median_s"], "(d/2 = 1.0)")
    print("ring Q", s1["median_Q"], "\nsquare Q", s2["median_Q"])
    assert s1["frac_degenerate"] == 0.0 and s2["frac_degenerate"] == 0.0
    assert abs(s1["median_Q"]["8"] - 0.5) <= 0.35 * 0.5 and abs(s2["median_Q"]["8"] - 1.0) <= 0.35 * 1.0
    assert abs(s1["median_s"]) <= 0.2 and abs(s2["median_s"]) <= 0.2


def test_clique_union_degenerate() -> None:
    a = np.zeros((100, 100))
    for b in range(10):
        a[b * 10 : (b + 1) * 10, b * 10 : (b + 1) * 10] = 1.0
    np.fill_diagonal(a, 0.0)
    _, s = _measure(a)  # giant component = una sola K10
    assert s["frac_degenerate"] > 0.5 or s["median_Q"]["8"] < 0.05


def test_random_regular_nondiffusive_values_reported() -> None:
    _, s = _measure(R.random_regular(500, 12, R.rng_from_key((7, 7, 7))))
    print("\nRR k12 median s", s["median_s"], "f_nd", s["f_nd"])
    assert s["median_s"] < 0.0 or s["frac_degenerate"] > 0.0
