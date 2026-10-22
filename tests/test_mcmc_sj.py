"""Tests de omega.causality.mcmc_sj (ES-0, docs/OMEGA_ES0_PILOTO.md §5)."""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import (
    init_general,
    init_two_order,
    measure,
    s_sj,
    sweep_general,
    sweep_two_order,
)
from omega.causality.orders import is_partial_order


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _n_rel(c: npt.NDArray[np.bool_]) -> tuple[float, int]:
    return float(c.sum()), 0


def _natural_orders(n: int) -> list[npt.NDArray[np.bool_]]:
    iu = np.triu_indices(n, k=1)
    out = []
    for bits in itertools.product([False, True], repeat=len(iu[0])):
        c = np.zeros((n, n), dtype=bool)
        c[iu] = bits
        cf = c.astype(np.int64)
        if not np.any(((cf @ cf) > 0) & ~c):
            out.append(c)
    return out


def _tv(emp: Counter[bytes], exact: dict[bytes, float]) -> float:
    tot = sum(emp.values())
    assert set(emp) <= set(exact)
    return 0.5 * sum(abs(emp.get(k, 0) / tot - p) for k, p in exact.items())


@pytest.mark.parametrize("lam", [0.7, -0.5])
def test_general_exact_distribution(lam: float) -> None:
    orders = _natural_orders(4)
    assert len(orders) == 40  # ordenes etiquetados naturalmente en 4 elementos
    w = {c.tobytes(): math.exp(lam * _n_rel(c)[0]) for c in orders}
    z = sum(w.values())
    exact = {k: x / z for k, x in w.items()}
    st = init_general(np.zeros((4, 4), dtype=bool), _n_rel)
    r = _rng(11 if lam > 0 else 12)
    emp: Counter[bytes] = Counter()
    for _ in range(40000):
        sweep_general(st, lam, r, _n_rel)
        emp[st.c.tobytes()] += 1
    assert _tv(emp, exact) < 0.03  # con lam desplazado en 0.2 el TV ronda 0.1


@pytest.mark.parametrize("lam", [0.7, -0.5])
def test_two_order_exact_distribution(lam: float) -> None:
    perms = [np.array(p, dtype=np.int64) for p in itertools.permutations(range(4))]
    w = {p.tobytes(): math.exp(lam * _n_rel(order_from_v(p))[0]) for p in perms}
    z = sum(w.values())
    exact = {k: x / z for k, x in w.items()}
    st = init_two_order(np.arange(4), _n_rel)
    r = _rng(21 if lam > 0 else 22)
    emp: Counter[bytes] = Counter()
    for _ in range(40000):
        sweep_two_order(st, lam, r, _n_rel)
        emp[st.v.tobytes()] += 1
    assert _tv(emp, exact) < 0.03  # con lam desplazado en 0.2 el TV ronda 0.1


def test_s_sj_antichain_and_random() -> None:
    assert s_sj(np.zeros((64, 64), dtype=bool)) == (0.0, 0)
    c = order_from_v(_rng(31).permutation(64).astype(np.int64))
    s, nr = s_sj(c)
    assert math.isfinite(s)
    assert 7 <= nr <= 9


def test_general_sampler_stays_ordered() -> None:
    c = order_from_v(_rng(41).permutation(16).astype(np.int64))
    st = init_general(c)
    r = _rng(42)
    for _ in range(3):
        cnt = sweep_general(st, 0.3, r)
        assert cnt.proposals == 120 and cnt.accepted <= cnt.valid <= cnt.proposals
        assert not np.any(np.tril(st.c))
        assert is_partial_order(st.c)
        assert st.s == pytest.approx(s_sj(st.c)[0])
    m = measure(st)
    assert 0.0 <= m["rel_frac"] <= 1.0 and m["height"] >= 1
