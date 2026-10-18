"""Omega-C0: dinamica de descenso proyectado."""

from __future__ import annotations

import numpy as np

from omega.c0.dynamics import evolve_c0
from omega.c0.functional import action_c0, params_from_targets

N = 30


def _start(seed: int) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(seed))
    u = np.triu(rng.random((N, N)) * 2 * 6 / (N - 1), k=1)
    return np.asarray(u + u.T, dtype=np.float64)


def test_monotone_decrease_box_and_symmetry() -> None:
    p = params_from_targets(2, 6, 1.0)
    w0 = _start(0)
    prev = action_c0(w0, p)
    for steps in (5, 50, 300):
        r = evolve_c0(w0, p, max_steps=steps)
        assert r["s_final"] <= prev + 1e-12 or steps == 5
        w = r["w"]
        assert r["s_final"] <= action_c0(w0, p)
        assert w.min() >= 0.0 and w.max() <= 1.0
        assert np.array_equal(w, w.T) and not w.diagonal().any()
        assert r["s_final"] == action_c0(w, p) or abs(r["s_final"] - action_c0(w, p)) < 1e-9
        prev = r["s_final"]


def test_equivariance_under_permutation() -> None:
    p = params_from_targets(0, 6, 1.0)
    w0 = _start(1)
    perm = np.random.Generator(np.random.PCG64(9)).permutation(N)
    r1 = evolve_c0(w0, p, max_steps=60)
    r2 = evolve_c0(w0[np.ix_(perm, perm)], p, max_steps=60)
    w1 = r1["w"][np.ix_(perm, perm)]
    assert np.max(np.abs(w1 - r2["w"])) <= 1e-8
