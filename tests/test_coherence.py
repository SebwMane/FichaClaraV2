"""Pruebas rapidas de omega/diagnostics/coherence.py y de los generadores nuevos de tools/coh_l0b.py."""

from __future__ import annotations

import numpy as np
from scipy.sparse.csgraph import connected_components

from omega.c0.references import rng_from_key
from omega.diagnostics.coherence import coherence_status, edge_coherence_profile
from tools import coh_l0b as C
from tools import p1d3_panel as P


def _status(adj: object, key: int = 1) -> tuple[str, dict[str, object]]:
    s = edge_coherence_profile(adj, rng_from_key((20261012, 999, key)))
    return coherence_status(s), s


def test_torus_square_coherente() -> None:
    st, s = _status(P.torus_lattice(60, 2))
    assert st == "COHERENTE", s
    assert 0.7 <= s["median_gamma"] <= 1.3  # type: ignore[operator]


def test_ring_coherente() -> None:
    st, s = _status(P.ring(2000, 2))
    assert st == "COHERENTE", s


def test_random_regular_not_coherente() -> None:
    st, s = _status(P.random_regular(3000, 6, rng_from_key((1, 2, 3))))
    assert st in ("INCOHERENTE", "SIN_VENTANA"), s


def test_caveman_degenerado() -> None:
    st, s = _status(P.caveman(2000, 12))
    assert st == "DEGENERADO", s


def test_status_table() -> None:
    base = {"r_w": 10, "median_gamma": 1.0, "f_deg": 0.0, "f_low": 0.0}
    assert coherence_status({**base, "r_w": 3}) == "SIN_VENTANA"
    assert coherence_status({**base, "f_deg": 0.6}) == "DEGENERADO"
    assert coherence_status(base) == "COHERENTE"
    assert coherence_status({**base, "median_gamma": 0.4}) == "INCOHERENTE"
    assert coherence_status({**base, "median_gamma": 0.6}) == "INTERMEDIO"


def test_heisenberg_degree_and_connected() -> None:
    a = C.heisenberg(9)
    assert a.shape == (729, 729)
    assert np.all(np.diff(a.indptr) == 4)
    assert connected_components(a, directed=False)[0] == 1


def test_new_generators() -> None:
    t = C.prufer_tree(500, rng_from_key((1, 1, 1)))
    assert t.nnz == 2 * 499 and connected_components(t, directed=False)[0] == 1
    s = C.subdivided_binary_tree(5, 8)
    assert s.shape[0] == 31 + 30 * 7 and s.nnz == 2 * 30 * 8
    assert connected_components(s, directed=False)[0] == 1
    c = C.husimi_cactus(1001)
    assert c.shape[0] == 1001 and c.nnz == 2 * 3 * 500
    assert connected_components(c, directed=False)[0] == 1
    p = C.cycle_times_rr(10, 50, 4, rng_from_key((1, 1, 2)))
    assert p.shape[0] == 500 and np.all(np.diff(p.indptr) == 6)
