"""Pruebas de config/seeds (D-28)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings import SeedConfig


def test_seed_key_structure_and_validation() -> None:
    sc = SeedConfig(master_entropy=99, experiment_id=3, replicates=4)
    assert seed_key(sc, 5, 2) == SeedKey(99, (3, 5, 2))
    with pytest.raises(ValueError):
        seed_key(sc, -1, 0)
    with pytest.raises(ValueError):
        seed_key(sc, 0, 4)
    with pytest.raises(ValueError):
        SeedKey(-1, ())


def test_determinism_bit_by_bit() -> None:
    k = SeedKey(42, (1, 2, 3))
    a = make_rng(k).random(1000)
    b = make_rng(SeedKey(42, (1, 2, 3))).random(1000)
    assert a.tobytes() == b.tobytes()


def test_different_keys_give_different_streams() -> None:
    sc = SeedConfig(master_entropy=1, experiment_id=0, replicates=5)
    draws = {
        (p, r): make_rng(seed_key(sc, p, r)).random(8).tobytes() for p in range(3) for r in range(5)
    }
    assert len(set(draws.values())) == len(draws)
    other = make_rng(SeedKey(2, (0, 0, 0))).random(8).tobytes()
    assert other not in draws.values()


def test_independent_of_creation_order() -> None:
    sc = SeedConfig(master_entropy=7, experiment_id=1)
    first = make_rng(seed_key(sc, 4, 1)).random(5)
    for p in range(10):
        make_rng(seed_key(sc, p, 0)).random(100)
    again = make_rng(seed_key(sc, 4, 1)).random(5)
    assert np.array_equal(first, again)


def test_generator_is_pcg64() -> None:
    g = make_rng(SeedKey(1, (0,)))
    assert isinstance(g, np.random.Generator)
    assert isinstance(g.bit_generator, np.random.PCG64)
