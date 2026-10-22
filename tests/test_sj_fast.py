"""Tests de omega.causality.sj_fast (R17-3, docs/OMEGA_CONSEJO_REV17_ES0.md §5.2): mismo valor que mcmc_sj.s_sj."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import (
    init_general,
    init_two_order,
    s_sj,
    sweep_general,
    sweep_two_order,
)
from omega.causality.orders import (
    kr_order,
    layered_2order,
    random_bipartite,
    transitive_percolation,
)
from omega.causality.sj_fast import (
    SJEvaluator,
    evaluate,
    s_sj_fast,
    sweep_general_fast,
    sweep_two_order_fast,
)
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

TOL = 1e-9
Bool = npt.NDArray[np.bool_]


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261102, spawn_key=(327, 99, seed))))


def _same(c: Bool) -> None:
    s0, nr0 = s_sj(c)
    s1, nr1 = s_sj_fast(c)
    assert nr1 == nr0
    assert abs(s1 - s0) <= TOL


def _kr3(n: int, rng: np.random.Generator) -> Bool:
    a = n // 4
    c = np.zeros((n, n), dtype=bool)
    c[:a, n - a :] = True
    c[:a, a : n - a] = rng.random((a, n - 2 * a)) < 0.5
    c[a : n - a, n - a :] = rng.random((n - 2 * a, a)) < 0.5
    return c


@pytest.mark.parametrize("n", [32, 64])
def test_stateless_matches_on_families(n: int) -> None:
    r = _rng(n)
    cases: list[Bool] = []
    for _ in range(4):
        cases.append(order_from_v(r.permutation(n).astype(np.int64)))
        cases.append(causal_matrix(sprinkle_diamond(n, 2, r)))
        cases.append(causal_matrix(sprinkle_diamond(n, 3, r)))
        cases.append(kr_order(n, r))
        cases.append(_kr3(n, r))
        cases.append(transitive_percolation(n, 0.15, r))
        cases.append(layered_2order(n, 4, r))
    cases.append(np.triu(np.ones((n, n), dtype=bool), k=1))  # cadena total
    has_region = 0
    for c in cases:
        _same(c)
        has_region += int(s_sj(c)[1] > 0)
    assert has_region >= len(cases) // 2


@pytest.mark.parametrize("n", [32, 64])
def test_region_empty_cases(n: int) -> None:
    r = _rng(100 + n)
    empty = [np.zeros((n, n), dtype=bool), random_bipartite(n, 0.5, r)]  # anticadena; altura 2 (|I| = 2)
    for c in empty:
        assert s_sj(c) == (0.0, 0)
        assert s_sj_fast(c) == (0.0, 0)
        res = evaluate(c)
        assert (res.s, res.n_r, res.fallback) == (0.0, 0, False)
        ev = SJEvaluator(c)
        assert (ev.s, ev.n_r) == (0.0, 0)


def test_fallback_is_bitwise_original() -> None:
    c = order_from_v(_rng(7).permutation(64).astype(np.int64))
    res = evaluate(c, guard=1e6)  # banda enorme: siempre se usa s_sj
    assert res.fallback
    assert (res.s, res.n_r) == s_sj(c)


def test_evaluator_moves_and_rollback_exact() -> None:
    n = 48
    r = _rng(8)
    ev = SJEvaluator.from_two_order(r.permutation(n).astype(np.int64))
    for _ in range(60):
        a, b = (int(x) for x in r.choice(n, size=2, replace=False))
        before = (ev.c.copy(), ev.m.copy(), ev.b.copy(), ev.s)
        res = ev.propose_swap(a, b)
        assert ev.check_state()
        s0, nr0 = s_sj(ev.c)
        assert res.n_r == nr0 and abs(res.s - s0) <= TOL
        if r.random() < 0.5:
            ev.rollback()
            assert np.array_equal(ev.c, before[0]) and np.array_equal(ev.m, before[1])
            assert np.array_equal(ev.b, before[2]) and ev.s == before[3]
        else:
            ev.commit()
    with pytest.raises(RuntimeError):
        ev.propose_swap(0, 1)
        ev.propose_swap(0, 1)
    ev.rollback()
    with pytest.raises(RuntimeError):
        ev.rollback()


def _check_two_order_chain(n: int, lam: float, seed: int, steps: int) -> int:
    """Cada propuesta (aceptada o no) comparada con s_sj; devuelve cuantas tuvieron region."""
    r = _rng(seed)
    ev = SJEvaluator.from_two_order(r.permutation(n).astype(np.int64))
    with_region = 0
    for _ in range(steps):
        a = int(r.integers(0, n))
        b = int(r.integers(0, n - 1))
        b += b >= a
        s_old = ev.s
        res = ev.propose_swap(a, b)
        s0, nr0 = s_sj(ev.c)
        assert res.n_r == nr0 and abs(res.s - s0) <= TOL
        with_region += int(nr0 > 0)
        x = lam * (res.s - s_old)
        if x >= 0 or r.random() < math.exp(x):
            ev.commit()
        else:
            ev.rollback()
    assert ev.check_state()
    return with_region


def _check_general_chain(c0: Bool, lam: float, seed: int, steps: int) -> tuple[int, int]:
    r = _rng(seed)
    ev = SJEvaluator(c0)
    c = ev.c
    n = c.shape[0]
    with_region = empty = 0
    for _ in range(steps):
        i, j = sorted(int(x) for x in r.choice(n, size=2, replace=False))
        if c[i, j]:
            if np.any(c[i] & c[:, j]):
                continue
        elif np.any(c[j] & ~c[i]) or np.any(c[:, i] & ~c[:, j]):
            continue
        s_old = ev.s
        res = ev.propose_toggle(i, j)
        s0, nr0 = s_sj(ev.c)
        assert res.n_r == nr0 and abs(res.s - s0) <= TOL
        if nr0 == 0:
            assert (res.s, res.n_r) == (0.0, 0)
            empty += 1
        else:
            with_region += 1
        x = lam * (res.s - s_old)
        if x >= 0 or r.random() < math.exp(x):
            ev.commit()
        else:
            ev.rollback()
    assert ev.check_state()
    return with_region, empty


@pytest.mark.parametrize(("n", "lam", "steps"), [(32, 0.0, 1000), (32, 28.0, 1000), (32, -28.0, 1000), (64, 7.0, 800), (64, -7.0, 800)])
def test_two_order_chain_every_proposal(n: int, lam: float, steps: int) -> None:
    assert _check_two_order_chain(n, lam, 200 + n + int(lam), steps) > steps // 2


@pytest.mark.parametrize("n", [32, 64])
def test_general_chain_every_proposal(n: int) -> None:
    r = _rng(300 + n)
    starts = [
        order_from_v(r.permutation(n).astype(np.int64)),
        _kr3(n, r),
        np.triu(np.ones((n, n), dtype=bool), k=1),
    ]
    tot_region = 0
    for k, c0 in enumerate(starts):
        lam = (0.0, 18.0, -18.0)[k]
        wr, _ = _check_general_chain(c0, lam, 310 + 10 * n + k, 1500 if n == 32 else 1000)
        tot_region += wr
    assert tot_region > 100
    # desde el orden vacio la region esta vacia al principio (convencion (0.0, 0))
    _, empty = _check_general_chain(np.zeros((n, n), dtype=bool), 0.0, 320 + n, 200)
    assert empty > 0


def test_sweeps_reproduce_mcmc_sj() -> None:
    """Mismas semillas: mismas trayectorias, conteos y S que los muestreadores de mcmc_sj."""
    n = 32
    v0 = _rng(400).permutation(n).astype(np.int64)
    st = init_two_order(v0)
    ev = SJEvaluator.from_two_order(v0)
    r1, r2 = _rng(401), _rng(401)
    for _ in range(2):
        assert sweep_two_order(st, -28.0, r1) == sweep_two_order_fast(ev, -28.0, r2)
        assert np.array_equal(st.v, ev.v) and abs(st.s - ev.s) <= TOL and st.n_r == ev.n_r
    c0 = _kr3(n, _rng(402))
    stg = init_general(c0)
    evg = SJEvaluator(c0)
    r1, r2 = _rng(403), _rng(403)
    for _ in range(2):
        assert sweep_general(stg, 17.0, r1) == sweep_general_fast(evg, 17.0, r2)
        assert np.array_equal(stg.c, evg.c) and abs(stg.s - evg.s) <= TOL and stg.n_r == evg.n_r
    assert evg.check_state() and ev.check_state()
