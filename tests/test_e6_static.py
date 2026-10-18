"""Tests de tools/e6_static.py: valores analiticos en Z^d, diamante, circulante, clases y rng."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import e6_static as e6  # noqa: E402
from omega.c0.references import rng_from_key  # noqa: E402
from p1d3_panel import torus_lattice  # noqa: E402


@pytest.mark.parametrize("d,side", [(1, 6), (2, 6), (3, 5), (4, 5)])
def test_torus_analytic(d: int, side: int) -> None:
    loc = e6.local_quantities(torus_lattice(side, d))
    assert (loc["deg"] == 2 * d).all()
    assert (loc["c4"] == 2 * d * (d - 1)).all()
    assert not loc["tri"].any()
    if d >= 2:
        assert np.allclose(loc["f4"], 2 * (d - 1) / (2 * d - 1))
    else:
        assert (loc["f4"] == 0).all()


def test_torus_edges_in_square() -> None:
    A = e6.local_quantities(torus_lattice(6, 2))["A"]
    u, v = e6.sample_edges(A, (1, 2, 3, 7))
    r = e6.edge_rest(A, u, v)
    assert r["C-SQ"].all() and not r["C-TRI"].any()
    A1 = e6.local_quantities(torus_lattice(8, 1))["A"]
    u, v = e6.sample_edges(A1, (1, 2, 3, 7))
    assert not e6.edge_rest(A1, u, v)["C-SQ"].any()
    assert u.size == 8


def test_diamond() -> None:
    A = e6.diamond_lattice(3)
    assert A.shape[0] == 8 * 27
    loc = e6.local_quantities(A)
    assert (loc["deg"] == 4).all()
    assert not loc["tri"].any()
    assert (loc["c4"] == 0).all()  # sin triangulos ni 4-ciclos => cintura >= 5; bipartito => 6
    A3 = (A @ A @ A).diagonal()
    assert (A3 == 0).all()
    A6 = np.linalg.matrix_power(A.toarray(), 6).diagonal()
    assert (A6 > 0).all()
    assert e6.diamond_lattice(10).shape[0] == 8000


def test_circulant_triangles() -> None:
    loc = e6.local_quantities(e6.circulant(50, (1, 2, 3)))
    assert (loc["deg"] == 6).all() and loc["tri"].all()


def test_lat_rule_on_torus() -> None:
    loc = e6.local_quantities(torus_lattice(5, 3))
    a = e6.rule_activity(loc, {"C-SQ": np.ones(3, bool), "C-TRI": np.zeros(3, bool)})
    assert a["P0-LAT(2)"] == 1.0 and a["P0-LAT(1)"] == 0.0 and a["P0-LAT(3)"] == 0.0
    assert a["P0-DEG(6)"] == 1.0 and a["P0-SQV(12)"] == 1.0 and a["P0-SQF(4/5)"] == 1.0


def prof(*v: float) -> dict[int, float]:
    return {d + 1: x for d, x in enumerate(v)}


def test_classes() -> None:
    c = e6.m_profile_classes
    assert c(prof(0.5, 0.55, 0.5, 0.5, 0.58, 0.5))["class"] == "CIEGO"
    assert c(prof(1, 1, 1, 1, 1, 1))["class"] == "CIEGO"
    assert c(prof(0, 0, 1, 0, 0, 0))["class"] == "SELECTOR"
    assert c(prof(0, 0.9, 0.2, 0, 0.5, 0))["class"] == "SELECTOR"
    assert c(prof(0, 1, 0, 1, 0, 0.3))["class"] == "SELECTOR"
    assert c(prof(1, 0, 0, 0, 0, 0))["class"] == "MONOTONO"
    assert c(prof(1, 1, 0.3, 0.1, 0, 0))["class"] == "MONOTONO"
    assert c(prof(0, 0, 0, 0, 0.5, 1))["class"] == "MONOTONO"
    assert c(prof(1, 0, 0, 0, 0, 1))["class"] == "INDETERMINADO"  # R no intervalo
    assert c(prof(1, 0.7, 0, 0, 0, 0))["class"] == "INDETERMINADO"  # d=2 fuera de R y de Q
    assert c(prof(0, 1, 1, 0, 0, 0))["class"] == "SELECTOR"
    assert c(prof(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))["class"] == "INDETERMINADO"  # R vacio
    # umbrales inclusivos
    assert c(prof(0.9, 0, 0, 0, 0, 0))["R"] == [1]
    assert c(prof(0.5, 0.5, 0.6, 0.5, 0.5, 0.5))["class"] == "CIEGO"  # max-min = 0.10 exacto
    assert c(prof(0.5, 0.5, 0.61, 0.5, 0.5, 0.5))["class"] != "CIEGO"


def test_dial() -> None:
    sel = e6.m_profile_classes(prof(0, 0, 1, 0, 0, 0))
    sel4 = e6.m_profile_classes(prof(0, 0, 0, 1, 0, 0))
    blind = e6.m_profile_classes(prof(1, 1, 1, 1, 1, 1))
    assert e6.dial_flag({"a": sel, "b": sel4})
    assert not e6.dial_flag({"a": sel, "b": sel})  # mismo R
    assert not e6.dial_flag({"a": blind, "b": e6.m_profile_classes(prof(0, 0, 0, 0, 0, 0))})
    mono = e6.m_profile_classes(prof(1, 0, 0, 0, 0, 0))
    mono2 = e6.m_profile_classes(prof(0, 0, 0, 0, 0, 1))
    assert not e6.dial_flag({"a": mono, "b": mono2})  # R distintos pero ninguno SELECTOR
    assert e6.dial_flag({"a": mono, "b": sel})


def test_rng_keys_and_ids() -> None:
    assert len(set(e6.GEOM_ID.values())) == len(e6.GEOM_ID)
    a = rng_from_key((e6.MASTER_E6, 13, 0)).random(5)
    b = rng_from_key((e6.MASTER_E6, 13, 0)).random(5)
    c = rng_from_key((e6.MASTER_E6, 13, 1)).random(5)
    assert (a == b).all() and not (a == c).all()
    A = e6.local_quantities(torus_lattice(40, 2))["A"]
    k = (e6.MASTER_E6, 10, 0, 7)
    assert all((x == y).all() for x, y in zip(e6.sample_edges(A, k), e6.sample_edges(A, k)))
    assert e6.sample_edges(A, k)[0].size == 3000
    assert set(e6._panel()) == set(e6.GEOM_ID)
