"""Tests del generador Kasner de vacio (R11-4): orden parcial, limite plano, densidad y K."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt
import pytest
from scipy import stats
from scipy.integrate import quad

from omega.causality.kasner import (
    in_kasner_interval,
    kasner_causal_matrix,
    kasner_directions,
    kasner_exponents,
    sprinkle_kasner_interval,
    support_table,
    validate_kasner_exponents,
)
from omega.causality.orders import is_partial_order
from omega.causality.sprinkling import causal_matrix

FLAT = (1.0, 0.0, 0.0)  # Kasner plano trivial: region de Milne de Minkowski
GENERIC = kasner_exponents(2.0)  # (-2/7, 3/7, 6/7), solo para tests


def _milne_to_minkowski(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    t, eta = x[:, 0], x[:, 1]
    return np.stack([t * np.cosh(eta), t * np.sinh(eta), x[:, 2], x[:, 3]], axis=1)


def test_vacuum_condition() -> None:
    for u in (1.0, 2.0, 7.5):
        p = np.array(kasner_exponents(u))
        assert abs(p.sum() - 1) < 1e-12 and abs((p * p).sum() - 1) < 1e-12
    validate_kasner_exponents(FLAT)
    with pytest.raises(ValueError):
        validate_kasner_exponents((1 / 3, 1 / 3, 1 / 3))  # sum p^2 = 1/3: no es vacio
    with pytest.raises(ValueError):
        validate_kasner_exponents((0.5, 0.5, 0.5))


def test_directions_symmetric_and_unit() -> None:
    d = kasner_directions(400)
    assert d.shape == (200, 3)
    assert np.allclose(np.linalg.norm(d, axis=1), 1.0)
    assert np.allclose(d[:3], np.eye(3))
    with pytest.raises(ValueError):
        kasner_directions(401)


def test_support_function_matches_quadrature() -> None:
    p = np.array(GENERIC)
    t = np.array([0.07, 0.3, 1.9])
    d = kasner_directions(20)
    h = support_table(t, p, d)
    for k in (0, 2, 5, 9):
        f = lambda s, u=d[k]: math.sqrt(float(np.sum(u * u * s ** (-2 * p))))  # noqa: E731
        ref = quad(f, t[0], t[2], epsabs=1e-13, epsrel=1e-13)[0]
        assert abs((h[2, k] - h[0, k]) - ref) < 1e-11
        assert np.all(np.diff(h[:, k]) > 0)


@pytest.mark.parametrize("p", [GENERIC, FLAT])
@pytest.mark.parametrize("k", [6, 400])
def test_partial_order(p: tuple[float, float, float], k: int) -> None:
    rng = np.random.default_rng(870101)
    x = sprinkle_kasner_interval(300, p, 0.3, 1.2, rng, k=k)
    c = kasner_causal_matrix(x, p, k)
    assert is_partial_order(c)
    assert 0.03 < c.mean() * 2 < 0.3  # fraccion de pares relacionados no trivial


def test_flat_limit_matches_minkowski() -> None:
    rng = np.random.default_rng(870102)
    x = sprinkle_kasner_interval(500, FLAT, 0.5, 1.5, rng, k=400)
    c = kasner_causal_matrix(x, FLAT, 400)
    cm = causal_matrix(_milne_to_minkowski(x))
    pairs = 500 * 499 / 2
    assert not np.any(cm & ~c)  # el poliedro contiene al cono verdadero: nunca faltan pares
    assert (c & ~cm).sum() / pairs < 2e-3  # medido ~5e-4 con K = 400
    assert cm.sum() / pairs > 0.05


def test_flat_limit_uniform_in_minkowski_diamond() -> None:
    """Con p = (1,0,0) el intervalo es el diamante de Minkowski entre T = t_lo y T = t_hi."""
    rng = np.random.default_rng(870103)
    t_lo, t_hi = 0.5, 1.5
    y = _milne_to_minkowski(sprinkle_kasner_interval(4000, FLAT, t_lo, t_hi, rng, k=400))
    h = 0.5 * (t_hi - t_lo)
    dt = y[:, 0] - t_lo
    assert np.all(np.linalg.norm(y[:, 1:], axis=1) <= np.minimum(dt, 2 * h - dt) * (1 + 0.01) + 1e-9)

    def cdf(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        z = np.clip(np.asarray(v, dtype=np.float64) / h, 0.0, 2.0)
        out: npt.NDArray[np.float64] = np.where(z <= 1, 0.5 * z**4, 1 - 0.5 * (2 - z) ** 4)
        return out

    assert stats.kstest(dt, cdf).pvalue > 1e-3


def test_density_proportional_to_t() -> None:
    """La densidad coordenada es prop. a t: frente a puntos uniformes en coordenadas del mismo
    intervalo, la ley de t del sprinkling es la uniforme reponderada por t."""
    p = GENERIC
    t_lo, t_hi = 0.2, 1.0
    rng = np.random.default_rng(870104)
    x = sprinkle_kasner_interval(6000, p, t_lo, t_hi, rng, k=200)
    assert np.all(in_kasner_interval(x, p, t_lo, t_hi, k=200))
    assert np.all((x[:, 0] >= t_lo) & (x[:, 0] <= t_hi))
    # referencia uniforme en coordenadas, mismo intervalo
    w = np.max(np.abs(x[:, 1:]), axis=0) * 1.05
    cand = np.column_stack([rng.uniform(t_lo, t_hi, 200000), rng.uniform(-1, 1, (200000, 3)) * w])
    ref = cand[in_kasner_interval(cand, p, t_lo, t_hi, k=200), 0]
    assert ref.size > 20000
    wt = ref / ref.sum()
    for power in (1, 2, -1):
        gr, gx = ref**power, x[:, 0] ** power
        mean_w = float(np.sum(wt * gr))
        se = float(np.std(gx) / math.sqrt(len(x)))
        assert abs(float(np.mean(gx)) - mean_w) < 4 * se
    # control de potencia: sin el peso t la media seria claramente menor
    se_t = float(np.std(x[:, 0]) / math.sqrt(len(x)))
    assert float(np.mean(x[:, 0])) - float(np.mean(ref)) > 10 * se_t


def test_convergence_in_k() -> None:
    rng = np.random.default_rng(870105)
    x = sprinkle_kasner_interval(400, GENERIC, 0.2, 1.0, rng, k=800)
    c200 = kasner_causal_matrix(x, GENERIC, 200)
    c800 = kasner_causal_matrix(x, GENERIC, 800)
    pairs = 400 * 399 / 2
    assert (c200 ^ c800).sum() / pairs < 5e-3  # medido ~1e-3
    assert c200.sum() >= c800.sum()  # mas direcciones: poliedro mas ajustado
