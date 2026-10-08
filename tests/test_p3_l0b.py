"""Tests rapidos de tools/p3_l0b.py (curvatura de Ollivier sobre grafos de referencia)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

import p3_l0b as P  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.landscape.references import connected_caveman  # noqa: E402


def test_torus_flat() -> None:
    r = P.curvature_sample(R.torus_lattice_3d(6), (1, 2, 3), max_edges=60)
    assert abs(r["median_kappa"]) <= 1e-9
    assert r["n_sampled"] == 60


def test_tree_negative() -> None:
    r = P.curvature_sample(R.random_tree(120, R.rng_from_key((1, 2, 3))), (1, 2, 3), max_edges=80)
    # con perezosidad 0.5 la mediana de un arbol es 0 (aristas de hoja: kappa = 0); el signo negativo esta en la cola
    assert r["median_kappa"] <= 1e-9
    assert r["f_neg"] > 0.3
    assert r["quantiles"]["10"] < -0.2


def test_caveman_positive() -> None:
    r = P.curvature_sample(connected_caveman(200, 8), (1, 2, 3), max_edges=80)
    assert r["median_kappa"] > 0
    assert r["f_pos"] > 0.5
