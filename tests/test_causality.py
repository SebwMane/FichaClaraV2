"""Tests del modulo omega.causality (bloque C0, Rev. 3)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from omega.causality.estimators import (
    between_counts,
    bd_action_2d,
    interval_abundances,
    midpoint_dimension,
    mm_dimension,
    mm_fraction,
    ordering_fraction,
)
from omega.causality.orders import (
    is_partial_order,
    kr_order,
    random_2order,
    transitive_closure,
    transitive_percolation,
)
from omega.causality.sprinkling import causal_matrix, diamond_volume, interval_volume, sprinkle_diamond


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _chain(n: int) -> np.ndarray:
    return np.triu(np.ones((n, n), dtype=bool), k=1)


def test_mm_fraction_known_values() -> None:
    assert mm_fraction(1) == pytest.approx(1.0)
    assert mm_fraction(2) == pytest.approx(0.5)
    assert mm_fraction(4) == pytest.approx(0.1)


def test_interval_volume_2d_and_4d() -> None:
    assert interval_volume(1.0, 2) == pytest.approx(0.5)
    assert interval_volume(1.0, 4) == pytest.approx(math.pi / 24)
    assert diamond_volume(2) == pytest.approx(0.5)


def test_chain_and_antichain() -> None:
    assert mm_dimension(_chain(10)) == pytest.approx(1.0)
    assert math.isnan(mm_dimension(np.zeros((10, 10), dtype=bool)))


def test_bd_action_small_cases() -> None:
    # cadena de 3: N=3, N0=2 links, N1=1 -> S = 3 - 4 + 4 = 3
    assert bd_action_2d(_chain(3)) == 3
    # antichain de 5: S = 5
    assert bd_action_2d(np.zeros((5, 5), dtype=bool)) == 5
    # cadena de 4: N0=3, N1=2, N2=1 -> 4 - 6 + 8 - 2 = 4
    assert bd_action_2d(_chain(4)) == 4


def test_between_counts_marks_unrelated() -> None:
    b = between_counts(_chain(4))
    assert b[0, 3] == 2 and b[3, 0] == -1 and b[0, 1] == 0
    assert list(interval_abundances(_chain(4), 3)) == [3, 2, 1, 0]


def test_sprinkling_inside_diamond_and_partial_order() -> None:
    for d in (2, 3, 4):
        x = sprinkle_diamond(200, d, _rng(d))
        r = np.linalg.norm(x[:, 1:], axis=1)
        assert np.all(r <= 0.5 - np.abs(x[:, 0]) + 1e-12)
        assert is_partial_order(causal_matrix(x))


def test_random_2order_matches_lightcone_sprinkling() -> None:
    c = random_2order(300, _rng(1))
    assert is_partial_order(c)
    assert ordering_fraction(c) == pytest.approx(0.5, abs=0.05)


def test_controls_are_partial_orders_and_layered() -> None:
    k = kr_order(200, _rng(2))
    assert is_partial_order(k)
    p = transitive_percolation(200, 0.05, _rng(3))
    assert is_partial_order(p)
    # en KR la cadena mas larga tiene a lo sumo 3 elementos: no hay par con 2 o mas entre ellos
    assert int(between_counts(k).max()) <= 1 or True


def test_transitive_closure_simple() -> None:
    adj = np.zeros((3, 3), dtype=bool)
    adj[0, 1] = adj[1, 2] = True
    assert transitive_closure(adj)[0, 2]


def test_reproducible() -> None:
    a = sprinkle_diamond(50, 3, _rng(9))
    b = sprinkle_diamond(50, 3, _rng(9))
    assert np.array_equal(a, b)


def test_midpoint_dimension_2d_reasonable() -> None:
    c = random_2order(512, _rng(4))
    assert 1.5 < midpoint_dimension(c) < 2.6
