"""Pruebas rapidas de tools/r3_0b_inspect.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import r3_0b as R  # noqa: E402
import r3_0b_inspect as I  # noqa: E402


def test_j3_prefix_consistency_and_monotone_L():
    key = (R.MASTER_R3, 8, 1, 0)
    big = R.poset_j3(30, R.rng_from_key(key))
    mask = 0
    prev = 0
    for n in (10, 15, 20, 25):
        small = R.poset_j3(n, R.rng_from_key(key))
        mask = (1 << n) - 1
        assert [b & mask for b in big[:n]] == small
        L = R.enumerate_downsets(R.make_poset("J3", 0, n, key), 10**9, False)["count"]
        assert L >= prev
        prev = L


def test_fit_synthetic():
    n = np.arange(5, 80)
    pts = [(int(k), int(round(np.exp(1.0 + 2.0 * np.sqrt(k))))) for k in n]
    g = I.growth_fits(pts)
    assert g["favors_sqrt"] and g["sqrt"]["sse"] < 1e-6 and abs(g["sqrt"]["b"] - 2.0) < 1e-3
    pts2 = [(int(k), int(round(k ** 3.0 * 2))) for k in n]
    g2 = I.growth_fits(pts2)
    assert not g2["favors_sqrt"] and abs(g2["log"]["b"] - 3.0) < 1e-3
    f = I.fit_sse(np.array([0., 1, 2]), np.array([1., 3, 5]))
    assert abs(f["a"] - 1) < 1e-9 and abs(f["b"] - 2) < 1e-9 and f["sse"] < 1e-12
