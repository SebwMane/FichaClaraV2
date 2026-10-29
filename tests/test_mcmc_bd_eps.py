"""Pruebas exactas de mcmc_bd_eps (ES-C0, OMEGA_ESC_PREREGISTRO.md §6): Delta S_eps incremental, B, balance
detallado exacto con N = 4 y equivalencia con mcmc_general para eps = 1.

Como en tests/test_bd_incremental.py, cada movimiento se ejecuta con el codigo real de ``sweep``: un generador de
prueba fuerza que la primera propuesta sea el par (i, j) con un u elegido, y que las demas se rechacen (beta = nan).
"""

from __future__ import annotations

import itertools
import math
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality import mcmc_general
from omega.causality.estimators import (
    bd_action_4d,
    bd_action_4d_smeared,
    between_counts,
)
from omega.causality.mcmc_bd_eps import action_eps, coef_table, init_state, sweep
from omega.causality.orders import is_partial_order, random_2order, transitive_closure
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

BoolMat = npt.NDArray[np.bool_]
EPS_SMEARED = (0.25, 0.6)


# ---------------------------------------------------------------- utilidades


class _Below:
    """Valor de u que hace verdadera la comparacion u < x (fuerza aceptar)."""

    def __lt__(self, other: object) -> bool:
        return True


class _U:
    """u[0] dado (o 'aceptar siempre'); el resto 0.5 (con beta = nan, u < exp(nan) es falso: se rechaza)."""

    def __init__(self, u0: float | None, rest: float) -> None:
        self.u0, self.rest = u0, rest

    def __getitem__(self, t: int) -> Any:
        if t == 0:
            return _Below() if self.u0 is None else self.u0
        return self.rest


class _OneMoveRng:
    """Imita los tres usos de rng en sweep: todas las propuestas son el par (i, j)."""

    def __init__(self, i: int, j: int, u0: float | None = None, rest: float = 0.5) -> None:
        self.i, self.j, self.u0, self.rest = i, j, u0, rest
        self.calls = 0

    def integers(self, lo: int, hi: int, size: int) -> npt.NDArray[np.int64]:
        self.calls += 1
        return np.full(size, self.i if self.calls == 1 else self.j - 1, dtype=np.int64)

    def random(self, size: int) -> _U:
        return _U(self.u0, self.rest)


def natural(c: BoolMat) -> BoolMat:
    perm = np.argsort(c.sum(axis=0), kind="stable")
    out: BoolMat = c[np.ix_(perm, perm)]
    assert not np.any(np.tril(out))
    return out


def full_between(c: BoolMat) -> npt.NDArray[np.int64]:
    cf = c.astype(np.int64)
    out: npt.NDArray[np.int64] = cf @ cf
    return out


def all_natural_orders(n: int) -> list[BoolMat]:
    pairs = list(itertools.combinations(range(n), 2))
    out = []
    for bits in itertools.product((False, True), repeat=len(pairs)):
        c = np.zeros((n, n), dtype=bool)
        for (a, b), v in zip(pairs, bits, strict=True):
            c[a, b] = v
        if np.array_equal(transitive_closure(c), c):
            out.append(c)
    return out


def expected_move(c: BoolMat, i: int, j: int) -> BoolMat | None:
    new = c.copy()
    if c[i, j]:
        if np.any(c[i] & c[:, j]):
            return None
        new[i, j] = False
    else:
        new[i, j] = True
        if not np.array_equal(transitive_closure(new), new):
            return None
    return new


def s_ref(c: BoolMat, eps: float) -> float:
    return float(bd_action_4d(c)) if eps == 1.0 else bd_action_4d_smeared(c, eps)


def forced_move(c: BoolMat, eps: float, i: int, j: int) -> tuple[bool, float, BoolMat, npt.NDArray[np.int64]]:
    """Ejecuta sweep con la propuesta (i, j) aceptada a la fuerza y las demas rechazadas (beta = nan).

    Devuelve (aceptada, Delta S incremental, C, B)."""
    st = init_state(c, eps)
    s0 = st.s
    acc, _ = sweep(st, float("nan"), _OneMoveRng(i, j))  # type: ignore[arg-type]
    return acc == 1, st.s - s0, st.c, st.b


def first_accepted(c: BoolMat, eps: float, i: int, j: int, u0: float, beta: float) -> bool:
    """Decision real de sweep sobre la primera propuesta (i, j) con u = u0 y beta finito.

    Las demas propuestas repiten el par con u = inf: un movimiento con a < 1 solo puede aceptarse en t = 0 (despues
    la inversa se acepta por el atajo y el directo vuelve a rechazarse). Asi, aceptada en t = 0 sii acc > 0."""
    st = init_state(c, eps)
    acc, _ = sweep(st, beta, _OneMoveRng(i, j, u0, float("inf")))  # type: ignore[arg-type]
    return acc > 0


# ---------------------------------------------------------------- pruebas


def test_table_matches_formula() -> None:
    g, const = coef_table(10, 1.0)
    assert const == 10.0 and list(g[:5]) == list(mcmc_general.COEF4) and not np.any(g[5:])
    for eps in EPS_SMEARED:
        g, const = coef_table(10, eps)
        k = np.arange(11, dtype=np.float64)
        s4 = (1 - eps) ** k * (1 - 9 * eps * k / (1 - eps) + 8 * eps**2 * k * (k - 1) / (1 - eps) ** 2
                               - 4 * eps**3 * k * (k - 1) * (k - 2) / (3 * (1 - eps) ** 3))
        np.testing.assert_allclose(g, -math.sqrt(eps) * eps * s4, rtol=1e-14, atol=1e-15)
        assert const == pytest.approx(math.sqrt(eps) * 10)
    with pytest.raises(ValueError):
        coef_table(4, 0.0)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_exhaustive_every_move(n: int) -> None:
    orders = all_natural_orders(n)
    assert len(orders) == {2: 2, 3: 7, 4: 40, 5: 357}[n]
    valid = 0
    for c in orders:
        for eps in (1.0, *EPS_SMEARED):
            assert init_state(c, eps).s == pytest.approx(s_ref(c, eps), abs=1e-12)
        for i, j in itertools.combinations(range(n), 2):
            exp = expected_move(c, i, j)
            for eps in (1.0, *EPS_SMEARED):
                ok, ds, cn, bn = forced_move(c, eps, i, j)
                assert ok == (exp is not None), (i, j, eps)
                if exp is None:
                    assert np.array_equal(cn, c) and ds == 0.0
                    continue
                assert np.array_equal(cn, exp)
                assert np.array_equal(bn, full_between(exp))
                full = s_ref(exp, eps) - s_ref(c, eps)
                if eps == 1.0:
                    assert ds == full and float(ds).is_integer()
                else:
                    assert abs(ds - full) <= 1e-12
            valid += exp is not None
    assert valid > 0


def _starts(rng: np.random.Generator) -> list[BoolMat]:
    return [
        natural(causal_matrix(sprinkle_diamond(16, 4, rng))),
        natural(random_2order(12, rng)),
        np.zeros((10, 10), dtype=bool),
        np.triu(np.ones((9, 9), dtype=bool), k=1),
    ]


@pytest.mark.parametrize("beta", [-1.0, -0.3, 0.0, 0.2, 1.0])
@pytest.mark.parametrize("eps", [1.0, 0.25, 0.6])
def test_chains_keep_action_and_between(eps: float, beta: float) -> None:
    rng = np.random.default_rng(int(1000 * (beta + 2)) + int(100 * eps))
    for c in _starts(rng):
        st = init_state(c, eps)
        accepted = 0
        for _ in range(5):
            acc, val = sweep(st, beta, rng)
            accepted += acc
            assert 0 <= acc <= val
            assert not np.any(np.tril(st.c))
            assert is_partial_order(st.c)
            assert np.array_equal(st.b, full_between(st.c))
            assert np.array_equal(st.b[st.c], between_counts(st.c)[st.c])
            if eps == 1.0:
                assert st.s == bd_action_4d(st.c)
            else:
                assert st.s == pytest.approx(bd_action_4d_smeared(st.c, eps), abs=1e-9)
            assert st.s == pytest.approx(action_eps(st.c, st.b, st.g, st.const), abs=1e-9)
        assert accepted > 0


@pytest.mark.parametrize("beta", [-0.5, 0.0, 0.05, 0.5])
def test_eps1_identical_to_mcmc_general(beta: float) -> None:
    rng0 = np.random.default_rng(7)
    for c in _starts(rng0):
        a = mcmc_general.init_state(c)
        e = init_state(c, 1.0)
        assert e.s == a.s
        ra, re = np.random.default_rng(99), np.random.default_rng(99)
        for _ in range(6):
            assert mcmc_general.sweep(a, beta, ra) == sweep(e, beta, re)
            assert np.array_equal(a.c, e.c) and np.array_equal(a.b, e.b) and a.s == e.s
        assert ra.random() == re.random()  # mismo flujo consumido


def _index4() -> tuple[list[BoolMat], dict[bytes, int]]:
    orders = all_natural_orders(4)
    return orders, {c.tobytes(): k for k, c in enumerate(orders)}


def _kernel(eps: float, beta: float) -> npt.NDArray[np.float64]:
    """Nucleo de una propuesta: par {i, j} con prob. 2/(N(N-1)) (i uniforme, j uniforme entre los otros)."""
    orders, index = _index4()
    n = 4
    p = np.zeros((len(orders), len(orders)))
    for x, c in enumerate(orders):
        for i, j in itertools.combinations(range(n), 2):
            ok, ds, cn, _ = forced_move(c, eps, i, j)
            if not ok:
                p[x, x] += 2 / (n * (n - 1))
                continue
            a = min(1.0, math.exp(-beta * ds))
            if a < 1.0:  # el sweep acepta justo cuando u < a
                assert first_accepted(c, eps, i, j, a * (1 - 1e-9), beta)
                assert not first_accepted(c, eps, i, j, a * (1 + 1e-9), beta)
            else:
                assert first_accepted(c, eps, i, j, 1 - 1e-12, beta)
            y = index[cn.tobytes()]
            p[x, y] += 2 / (n * (n - 1)) * a
            p[x, x] += 2 / (n * (n - 1)) * (1 - a)
    return p


def _pi(eps: float, beta: float) -> npt.NDArray[np.float64]:
    orders, _ = _index4()
    s = np.array([s_ref(c, eps) for c in orders])
    w = np.exp(-beta * (s - s.min()))
    out: npt.NDArray[np.float64] = w / w.sum()
    return out


@pytest.mark.parametrize(("eps", "beta"), [(1.0, -0.5), (1.0, 0.7), (0.25, -3.0), (0.25, 2.0), (0.6, -1.5)])
def test_exact_detailed_balance_n4(eps: float, beta: float) -> None:
    p = _kernel(eps, beta)
    pi = _pi(eps, beta)
    np.testing.assert_allclose(p.sum(axis=1), 1.0, atol=1e-14)
    flow = pi[:, None] * p
    np.testing.assert_allclose(flow, flow.T, atol=1e-15)
    np.testing.assert_allclose(pi @ p, pi, atol=1e-15)
    # irreducible: el nucleo conecta los 40 ordenes
    reach = (p > 0).astype(np.float64)
    assert np.all(np.linalg.matrix_power(reach + np.eye(40), 40) > 0)


def test_long_chain_tv_n4_eps025() -> None:
    eps, beta = 0.25, -3.0
    pi = _pi(eps, beta)
    assert 0.5 * np.abs(pi - 1 / 40).sum() > 0.15  # peso visiblemente no uniforme
    _, index = _index4()
    st = init_state(np.zeros((4, 4), dtype=bool), eps)
    rng = np.random.default_rng(20261019)
    counts = np.zeros(40)
    for _ in range(60_000):
        sweep(st, beta, rng)
        counts[index[st.c.tobytes()]] += 1
    tv = 0.5 * float(np.abs(counts / counts.sum() - pi).sum())
    assert tv < 0.02, tv
    assert st.s == pytest.approx(bd_action_4d_smeared(st.c, eps), abs=1e-9)
