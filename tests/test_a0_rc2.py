"""Pruebas rapidas de tools/a0_rc2.py (generadores nuevos y condiciones RC-2)."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import a0_rc2 as A  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402


def _deg(a: Any) -> np.ndarray:
    return np.asarray(a.sum(axis=1)).ravel()


def _nedges(a: Any) -> int:
    return int(a.nnz // 2)


def test_honeycomb_degree3_girth6() -> None:
    a = A.honeycomb_torus(20, 22)
    assert (_deg(a) == 3).all()
    assert A.locality_fraction(a, rng_from_key((1, 2, 3)), n_edges=10**9) == 0.0
    assert (a @ a @ a).diagonal().sum() == 0  # sin triangulos


def test_triangular_torus_degree6() -> None:
    a = A.triangular_torus(15)
    assert a.shape[0] == 225 and (_deg(a) == 6).all()


def test_apollonian_edges() -> None:
    n = 500
    a = A.apollonian(n, rng_from_key((1, 2, 3)))
    assert a.shape[0] == n and _nedges(a) == 3 * n - 6


def test_three_tree_edges() -> None:
    n = 500
    a = A.random_three_tree(n, rng_from_key((1, 2, 4)))
    assert a.shape[0] == n and _nedges(a) == 3 * n - 6


def test_other_generators_basic() -> None:
    assert (_deg(A.tree_of_grids(7, 10)) > 0).all() and A.tree_of_grids(7, 10).shape[0] == 700
    ba = A.barabasi_albert(300, 6, rng_from_key((1, 2, 5)))
    assert _nedges(ba) == 21 + 6 * (300 - 7)
    t = A.torus_with_diagonals(7, 0.2, rng_from_key((1, 2, 6)))
    assert _nedges(t) == 3 * 343 + round(0.2 * 3 * 343)
    assert A.lattice_open(5, 3).nnz // 2 == 3 * 4 * 25
    c = A.cartesian(A.cycle(10), A.path(4))
    assert c.shape[0] == 40 and _nedges(c) == 10 * 3 + 4 * 10
    cyl = A.rgg_cylinder(800, 12.0, rng_from_key((1, 2, 7)))
    assert cyl.shape[0] == 800


def _row(**kw: Any) -> dict[str, Any]:
    r: dict[str, Any] = {"locality": 0.9, "annulus_status": "CONEXO", "annulus": {"scales": [2, 4, 8]},
                         "coherence_status": "COHERENTE", "curvature": {"median": 0.1, "f_neg": 0.1}}
    r.update(kw)
    return r


def test_rc2_conditions() -> None:
    assert A.rc2_conditions(_row())["nondegenerate"]
    assert not A.rc2_conditions(_row(locality=0.3))["i_locality"]
    assert not A.rc2_conditions(_row(annulus_status="RAMIFICADO"))["ii_one_end"]
    assert not A.rc2_conditions(_row(annulus={"scales": [4]}))["ii_one_end"]
    assert A.rc2_conditions(_row(coherence_status="INTERMEDIO"))["iii_coherence"]
    assert not A.rc2_conditions(_row(coherence_status="ROTA"))["iii_coherence"]
    assert not A.rc2_conditions(_row(curvature={"median": 0.5, "f_neg": 0.1}))["iv_curvature"]
    assert not A.rc2_conditions(_row(curvature={"median": 0.1, "f_neg": 0.5}))["iv_curvature"]
    assert not A.rc2_conditions(_row(curvature={"median": 0.1, "f_neg": 0.5}))["nondegenerate"]
