"""Tests de tools/e6_static_v11.py: operadores de transporte, S2, C-5CIC, muestreo."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import e6_static as e6  # noqa: E402
import e6_static_v11 as v11  # noqa: E402
from p1d3_panel import _from_edges, torus_lattice  # noqa: E402


def adjsets(A):
    return [set(A.indices[A.indptr[i]:A.indptr[i + 1]].tolist()) for i in range(A.shape[0])]


def c4_tri(A):
    loc = e6.local_quantities(A)
    return loc


@pytest.mark.parametrize("d,side", [(2, 6), (3, 5), (4, 5)])
def test_truncation(d: int, side: int) -> None:
    A = v11.truncation(side, d)
    assert A.shape[0] == 2 * d * side**d
    loc = e6.local_quantities(A)
    assert (loc["deg"] == 3).all() and not loc["tri"].any()
    if d >= 3:
        assert (loc["c4"] == 0).all()
    else:
        assert (loc["c4"] > 0).all()  # gadget C4 en d=2


@pytest.mark.parametrize("d,side", [(1, 8), (2, 6), (3, 5), (4, 5)])
def test_line_graph(d: int, side: int) -> None:
    A = v11.line_graph(side, d)
    assert A.shape[0] == d * side**d
    assert (np.diff(A.indptr) == 4 * d - 2).all()
    assert (A != A.T).nnz == 0


@pytest.mark.parametrize("d,side", [(1, 6), (2, 6), (3, 5), (4, 5)])
def test_prod_k2(d: int, side: int) -> None:
    A = v11.prod_k2(side, d)
    assert A.shape[0] == 2 * side**d
    loc = e6.local_quantities(A)
    assert (loc["deg"] == 2 * d + 1).all()
    assert (loc["c4"] == 2 * d * d).all()


@pytest.mark.parametrize("d", [1, 2, 3, 4])
def test_s2(d: int) -> None:
    allv = lambda A: np.arange(A.shape[0])  # noqa: E731
    A = torus_lattice(7, d)
    assert (v11.s2_sizes(A, allv(A)) == 2 * d * d).all()
    B = v11.prod_k2(7, d)
    assert (v11.s2_sizes(B, allv(B)) == 2 * d * d + 2 * d).all()


def _cyc(n):
    i = np.arange(n)
    return _from_edges(n, i, (i + 1) % n)


def test_5cic() -> None:
    def frac(A):
        adj = adjsets(A)
        return [v11.in_5cycle(adj, u, v) for u in range(A.shape[0]) for v in adj[u] if u < v]
    assert all(frac(_cyc(5)))
    assert not any(frac(_cyc(6)))
    k4 = _from_edges(4, *map(np.array, zip(*[(a, b) for a in range(4) for b in range(a + 1, 4)])))
    assert not any(frac(k4))
    k5 = _from_edges(5, *map(np.array, zip(*[(a, b) for a in range(5) for b in range(a + 1, 5)])))
    assert all(frac(k5))
    # C5 con una cuerda: todas las aristas del ciclo siguen en el 5-ciclo, la cuerda no
    A = _from_edges(5, np.array([0, 1, 2, 3, 4, 0]), np.array([1, 2, 3, 4, 0, 2]))
    adj = adjsets(A)
    assert v11.in_5cycle(adj, 0, 1) and v11.in_5cycle(adj, 3, 4)
    assert not v11.in_5cycle(adj, 0, 2)  # el unico ciclo hamiltoniano es C5; la cuerda 0-2 no esta en ningun 5-ciclo


def test_subset_matches_full() -> None:
    A = e6.diamond_lattice(3)
    full = e6.local_quantities(A)
    verts = np.array([0, 5, 17, 100])
    sub = v11.local_quantities_subset(A, verts)
    for k in ("deg", "tri", "c4", "f4"):
        assert np.allclose(full[k][verts], sub[k])
    A = e6.local_quantities(e6.circulant(40, (1, 2, 3)))["A"]
    full = e6.local_quantities(A)
    sub = v11.local_quantities_subset(A, np.arange(40))
    for k in ("deg", "tri", "c4", "f4"):
        assert np.allclose(full[k], sub[k])


def test_ids_and_sampling_determinism() -> None:
    assert len(set(v11.GEOM_ID.values())) == len(v11.GEOM_ID)
    assert all(v11.GEOM_ID[k] == v for k, v in e6.GEOM_ID.items())
    assert set(v11.panel_v11()) == set(v11.GEOM_ID)
    from omega.c0.references import rng_from_key
    k = (e6.MASTER_E6, v11.GEOM_ID["L_Z6_5"], 0, 8)
    a = np.sort(rng_from_key(k).choice(93750, size=3000, replace=False))
    b = np.sort(rng_from_key(k).choice(93750, size=3000, replace=False))
    assert (a == b).all() and np.unique(a).size == 3000


def test_run_small_geometry_deterministic() -> None:
    P = v11.panel_v11()
    dim, ctor, rand = P["Z2K2_70"]
    g1 = v11.run_geometry("Z2K2_70", dim, ctor, rand)
    g2 = v11.run_geometry("Z2K2_70", dim, ctor, rand)
    assert g1["a"] == g2["a"]
    assert g1["a"]["P0-S2(8)"] == 0.0
