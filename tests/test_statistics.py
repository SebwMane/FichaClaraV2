"""Pruebas de statistics/summary."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.seeds import SeedKey, make_rng
from omega.statistics.summary import bootstrap_ci, proportion_ci


def test_bootstrap_ci_contains_mean_and_is_deterministic() -> None:
    v = make_rng(SeedKey(1, (0,))).normal(5.0, 1.0, size=200)
    lo, hi = bootstrap_ci(v, make_rng(SeedKey(2, (0,))))
    lo2, hi2 = bootstrap_ci(v, make_rng(SeedKey(2, (0,))))
    assert (lo, hi) == (lo2, hi2)
    assert lo < v.mean() < hi and hi - lo < 0.6


def test_bootstrap_ci_constant_and_validation(rng: np.random.Generator) -> None:
    assert bootstrap_ci(np.full(10, 3.0), rng) == (3.0, 3.0)
    with pytest.raises(ValueError):
        bootstrap_ci(np.array([]), rng)
    with pytest.raises(ValueError):
        bootstrap_ci(np.array([1.0, np.nan]), rng)
    with pytest.raises(ValueError):
        bootstrap_ci(np.ones(3), rng, level=1.0)


def test_input_not_mutated(rng: np.random.Generator) -> None:
    v = np.arange(10.0)
    before = v.copy()
    bootstrap_ci(v, rng)
    assert np.array_equal(v, before)


def test_wilson_known_values() -> None:
    lo, hi = proportion_ci(5, 10)
    assert lo == pytest.approx(0.2366, abs=1e-3) and hi == pytest.approx(0.7634, abs=1e-3)
    lo0, hi0 = proportion_ci(0, 10)
    assert lo0 == pytest.approx(0.0, abs=1e-12) and 0.25 < hi0 < 0.32
    lo1, hi1 = proportion_ci(10, 10)
    assert hi1 == pytest.approx(1.0, abs=1e-12) and lo1 > 0.68


def test_wilson_validation() -> None:
    with pytest.raises(ValueError):
        proportion_ci(11, 10)
    with pytest.raises(ValueError):
        proportion_ci(0, 0)
    with pytest.raises(TypeError):
        proportion_ci(1.0, 10)  # type: ignore[arg-type]
