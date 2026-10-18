"""Omega-C0: criterio de localidad y clasificador C0-L4."""

from __future__ import annotations

import numpy as np

from omega.c0.locality import classify_c0, degree_preserving_rewire, h_null, short_cycle_fraction, validate_locality
from omega.c0.references import random_tree, torus_lattice_3d


def _rng(s: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(s))


def test_short_cycle_fraction_known_cases() -> None:
    assert short_cycle_fraction(torus_lattice_3d(6)) == 1.0
    assert short_cycle_fraction(random_tree(60, _rng(1))) == 0.0
    tri = np.ones((3, 3)) - np.eye(3)
    assert short_cycle_fraction(tri) == 1.0
    c5 = np.zeros((5, 5))
    for i in range(5):
        c5[i, (i + 1) % 5] = c5[(i + 1) % 5, i] = 1.0
    assert short_cycle_fraction(c5) == 0.0


def test_rewire_preserves_degrees_and_is_deterministic() -> None:
    a = torus_lattice_3d(4) > 0
    r1 = degree_preserving_rewire(a, 500, _rng(5))
    r2 = degree_preserving_rewire(a, 500, _rng(5))
    assert np.array_equal(r1, r2)
    assert np.array_equal(r1.sum(axis=1), a.sum(axis=1))
    assert np.array_equal(r1, r1.T) and not r1.diagonal().any()


def test_h_null_deterministic() -> None:
    a = torus_lattice_3d(5)
    assert h_null(a, _rng(2)) == h_null(a, _rng(2))
    assert h_null(a, _rng(2)) > 1.15


def test_classify_empty_and_dense() -> None:
    assert classify_c0(np.zeros((20, 20)), _rng(0))["class"] == "VACIO"
    u = (np.ones((20, 20)) - np.eye(20)) * 0.3
    assert classify_c0(u, _rng(0))["class"] == "DENSO_TRIVIAL"


def test_validate_locality_ok() -> None:
    v = validate_locality(0)
    assert v["ok"], {k: r["class"] for k, r in v["graphs"].items()}
    assert len(v["graphs"]) == 9
