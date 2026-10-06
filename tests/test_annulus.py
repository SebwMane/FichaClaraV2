"""Pruebas rapidas de omega/diagnostics/annulus.py (docs/OMEGA_CIC_L0.md §3)."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from omega.c0.references import rng_from_key
from omega.diagnostics.annulus import annulus_profile, annulus_status
from tools import coh_l0b as C
from tools import p1d3_panel as P


def _st(adj: object, key: int = 1) -> tuple[str, dict[str, object]]:
    s = annulus_profile(adj, rng_from_key((20261013, 999, key)))
    return annulus_status(s), s


def test_torus_conexo() -> None:
    st, s = _st(P.torus_lattice(40, 2))
    assert st == "CONEXO", s


def test_ring_dos_extremos() -> None:
    st, s = _st(P.ring(2000, 4))
    assert st == "DOS_EXTREMOS", s


def test_trees_ramificado() -> None:
    st, s = _st(C.subdivided_binary_tree(8, 4))
    assert st == "RAMIFICADO", s


def test_complete_binary_tree_ramificado() -> None:
    n = 2**12 - 1
    child = np.arange(1, n, dtype=np.int64)
    st, s = _st(P._from_edges(n, (child - 1) // 2, child))
    assert st == "RAMIFICADO", s


def test_sin_ventana() -> None:
    k = 60
    iu, ju = np.triu_indices(k, 1)
    assert _st(P._from_edges(k, iu.astype(np.int64), ju.astype(np.int64)))[0] == "SIN_VENTANA"
    assert _st(P.random_regular(2000, 12, rng_from_key((1, 2, 3))))[0] == "SIN_VENTANA"


def test_status_table() -> None:
    base = {"f_med": 0.5, "c_med": 2.0, "f2_med": 0.97}
    assert annulus_status({**base, "f_med": None}) == "SIN_VENTANA"
    assert annulus_status({**base, "f_med": 0.95}) == "CONEXO"
    assert annulus_status(base) == "DOS_EXTREMOS"
    assert annulus_status({**base, "c_med": 3.0}) == "RAMIFICADO"
    assert isinstance(sparse.csr_array((2, 2)), sparse.csr_array)
