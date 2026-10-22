"""Auditoria Rev. 17 §8 paso 2: accion BD 4D completa frente a la actualizacion incremental de mcmc_general.

Cada movimiento se ejecuta con el codigo real de ``mcmc_general.sweep``: un generador de prueba fuerza que la
primera propuesta sea el par (i, j) y se acepte, y que las demas se rechacen. Asi se compara el Delta S
incremental de ``sweep`` con bd_action_4d(despues) - bd_action_4d(antes), exacto en enteros.

No existe muestreador incremental de la accion suavizada (bd_action_4d_smeared); aqui no se prueba.
"""

from __future__ import annotations

import itertools
import zlib
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.estimators import bd_action_4d
from omega.causality.mcmc_general import COEF4, action4, init_state, sweep
from omega.causality.orders import (
    is_partial_order,
    random_2order,
    transitive_closure,
    transitive_percolation,
)
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond
from omega.causality.torus import retarded_action_4d

BoolMat = npt.NDArray[np.bool_]


# ---------------------------------------------------------------- utilidades


class _Accept:
    """Valor de u que hace verdadera la comparacion u < exp(...) (fuerza la aceptacion)."""

    def __lt__(self, other: object) -> bool:
        return True


class _U:
    """u[0] fuerza aceptar; el resto es 0.5 (con beta = nan, u < exp(nan) es falso: se rechaza)."""

    def __getitem__(self, t: int) -> Any:
        return _Accept() if t == 0 else 0.5


class _OneMoveRng:
    """Imita los tres usos de rng en sweep: todas las propuestas son el par (i, j)."""

    def __init__(self, i: int, j: int) -> None:
        self.i, self.j = i, j
        self.calls = 0

    def integers(self, lo: int, hi: int, size: int) -> npt.NDArray[np.int64]:
        self.calls += 1
        # sweep: i_arr = integers(0, n); j_arr = integers(0, n - 1) + (j_arr >= i_arr)
        return np.full(size, self.i if self.calls == 1 else self.j - 1, dtype=np.int64)

    def random(self, size: int) -> _U:
        return _U()


def natural(c: BoolMat) -> BoolMat:
    """Reetiqueta con una extension lineal: si a < b, a tiene estrictamente menos ancestros que b."""
    perm = np.argsort(c.sum(axis=0), kind="stable")
    out: BoolMat = c[np.ix_(perm, perm)]
    assert not np.any(np.tril(out))
    return out


def full_between(c: BoolMat) -> npt.NDArray[np.int64]:
    cf = c.astype(np.int64)
    out: npt.NDArray[np.int64] = cf @ cf
    return out


def expected_move(c: BoolMat, i: int, j: int) -> BoolMat | None:
    """Movimiento de P4 (Rev. 3 §21) definido de forma independiente: None si no es valido."""
    new = c.copy()
    if c[i, j]:
        if np.any(c[i] & c[:, j]):
            return None  # relacion que no es link
        new[i, j] = False
    else:
        new[i, j] = True
        if not np.array_equal(transitive_closure(new), new):
            return None  # agregarla exige mas relaciones
    assert is_partial_order(new)
    return new


def check_move(c: BoolMat, i: int, j: int) -> bool:
    """Ejecuta sweep con un solo movimiento forzado y verifica Delta S y B. Devuelve si fue valido."""
    st = init_state(c)
    s0 = bd_action_4d(c)
    assert st.s == s0
    acc, _ = sweep(st, float("nan"), _OneMoveRng(i, j))  # type: ignore[arg-type]
    exp = expected_move(c, i, j)
    assert acc == (exp is not None), (i, j)
    if exp is None:
        assert np.array_equal(st.c, c) and st.s == s0
        return False
    assert np.array_equal(st.c, exp)
    ds_full = bd_action_4d(exp) - s0
    assert st.s - s0 == ds_full  # enteros exactos en coma flotante
    assert float(st.s).is_integer()
    assert np.array_equal(
        st.b, full_between(exp)
    )  # B mantenida en todos los pares, no solo los relacionados
    return True


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


def random_orders() -> list[tuple[str, BoolMat]]:
    rng = np.random.default_rng(20261019)
    out: list[tuple[str, BoolMat]] = []
    for n in (8, 12, 20, 30):
        for d in (2, 4):
            x = sprinkle_diamond(n, d, rng)
            out.append((f"sprinkling{d}d_N{n}", natural(causal_matrix(x))))
        out.append((f"2orden_N{n}", natural(random_2order(n, rng))))
        for p in (0.1, 0.3):
            out.append(
                (f"percolacion{p}_N{n}", natural(transitive_percolation(n, p, rng)))
            )
        up = np.triu(rng.random((n, n)) < 0.15, k=1)
        out.append((f"transitivo_N{n}", transitive_closure(up)))
    return out


RANDOM = random_orders()


# ---------------------------------------------------------------- pruebas


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_exhaustive_small_orders_every_move(n: int) -> None:
    orders = all_natural_orders(n)
    assert (
        len(orders) == {2: 2, 3: 7, 4: 40, 5: 357}[n]
    )  # 40 coincide con P4-0 (Rev. 3 §21)
    valid = 0
    for c in orders:
        for i, j in itertools.combinations(range(n), 2):
            valid += check_move(c, i, j)
    assert valid > 0


@pytest.mark.parametrize(("name", "c"), RANDOM, ids=[r[0] for r in RANDOM])
def test_random_orders_moves(name: str, c: BoolMat) -> None:
    n = c.shape[0]
    rng = np.random.default_rng(zlib.crc32(name.encode()))
    pairs = list(itertools.combinations(range(n), 2))
    links = [p for p in pairs if c[p] and not np.any(c[p[0]] & c[:, p[1]])]
    addable = [p for p in pairs if not c[p] and expected_move(c, *p) is not None]
    taken = set(links) | set(addable)
    other = [p for p in pairs if p not in taken]
    picks: list[tuple[int, int]] = []
    for group, k in ((links, 12), (addable, 12), (other, 6)):
        if group:
            idx = rng.choice(len(group), size=min(k, len(group)), replace=False)
            picks += [group[t] for t in idx]
    assert links  # todo orden no vacio tiene links
    for i, j in picks:
        check_move(c, i, j)


@pytest.mark.parametrize(("name", "c"), RANDOM, ids=[r[0] for r in RANDOM])
def test_three_full_actions_agree(name: str, c: BoolMat) -> None:
    b = full_between(c)
    s = bd_action_4d(c)
    assert action4(c, b) == s
    assert retarded_action_4d(c, 1.0).sum() == pytest.approx(s, abs=1e-9)
    # COEF4[k] es el aporte de un par con k elementos entre; k >= 4 aporta 0
    assert list(COEF4) == [-1.0, 9.0, -16.0, 8.0, 0.0]


@pytest.mark.parametrize("beta", [-0.5, -0.1, 0.0, 0.05, 0.5])
def test_sweep_keeps_action_and_between_exact(beta: float) -> None:
    rng = np.random.default_rng(int(1000 * (beta + 1)))
    starts = [
        natural(causal_matrix(sprinkle_diamond(16, 4, rng))),
        natural(random_2order(12, rng)),
        np.zeros((10, 10), dtype=bool),
    ]
    for c in starts:
        st = init_state(c)
        accepted = 0
        for _ in range(4):
            acc, _ = sweep(st, beta, rng)
            accepted += acc
            assert not np.any(np.tril(st.c))
            assert is_partial_order(st.c)
            assert st.s == bd_action_4d(st.c)
            assert np.array_equal(st.b, full_between(st.c))
        assert accepted > 0
