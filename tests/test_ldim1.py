"""Tests de L-DIM-1: caracteristicas ciegas y logica del evaluador (docs/OMEGA_LDIM1_PRERREGISTRO.md)."""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import ldim1_eval as ev  # noqa: E402
import ldim1_features as lf  # noqa: E402
from p1d3_panel import ring  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402


def test_clique_phi_zero_psi_inverse():
    n = 200
    f = lf.extract(lf.clique(n), "x", rng_from_key((1, 2, 3)))
    assert f["phi"] == 0.0
    assert abs(f["psi"] - 1.0 / (n - 1)) < 1e-12


def test_cycle_phi():
    f = lf.extract(ring(4000, 2), "x", rng_from_key((1, 2, 3)))
    r = f["r_star"]
    assert r > 1
    assert abs(f["phi"] - 2.0 / (2 * r + 1)) < 1e-9
    assert abs(f["psi"] - 0.5) < 1e-12


def test_blindness_signature_and_id():
    assert list(inspect.signature(lf.extract).parameters) == ["adj", "gid", "rng"]
    key = (lf.MASTER, 1, 0)
    assert len(lf.graph_id(key)) == 12 and lf.graph_id(key) == lf.graph_id(key)
    assert set(lf.extract(ring(300, 2), "abc", rng_from_key(key))) == {"id", "r_star", "max_window", "n_giant", "m", "phi", "psi"}
    ids = [s["id"] for s in lf.specs(False)]
    assert len(set(ids)) == len(ids)


def test_run_length_logic():
    ok = np.array([0, 1, 1, 1, 1, 1, 0, 1, 1], dtype=bool)
    dst = np.array([0, 3, 3, 3, 3, 3, 0, 3, 3])
    assert ev.longest_run(ok, dst) == (1, 5, 3)
    assert ev.l3_pass(ok, dst)  # 5 puntos = 1 decada
    assert not ev.l3_pass(np.array([1, 1, 1, 1], dtype=bool), np.array([2, 2, 2, 2]))  # 4 puntos = 0.75
    assert not ev.l3_pass(np.array([1, 1, 1, 1, 1, 1], dtype=bool), np.array([2, 2, 3, 3, 3, 3]))  # cambio de d*
    assert ev.longest_run(np.zeros(5, dtype=bool), np.zeros(5, dtype=int)) == (0, 0, 0)


def test_dstar_and_l1_l2_synthetic():
    d = np.array([1, 2, 2, 3, 3, 4])
    st = np.array(["reticulo", "reticulo", "k8", "k8", "k16", "k16"])
    F = np.array([[5.0, 2.0, 2.0, 1.0, 1.0, 3.0], [1.0, 2.0, 2.0, 3.0, 3.0, 4.0]])
    pooled, l1, l2 = ev.l1_l2(F, d, st)
    assert pooled.tolist() == [3, 1]
    assert l1.tolist() == [True, False]
    assert l2[0] is np.False_ or not l2[0]  # estrato reticulo solo contiene d = 1, 2 -> d* = 2 != 3
