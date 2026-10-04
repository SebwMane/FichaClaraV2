"""WP-D: diagnosticos de ensemble sobre series sinteticas con tau_int conocido."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.settings11 import EnsembleConfig
from omega.contracts import ChainResult
from omega.config.settings11 import Engine
from omega.statistics.ensemble import (
    autocorrelation,
    effective_sample_size,
    ensemble_summary,
    geweke_z,
    integrated_autocorr_time,
    post_burn_in,
    split_rhat,
)


def ar1(n: int, phi: float, seed: int, mean: float = 0.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    e = rng.standard_normal(n)
    x = np.empty(n)
    x[0] = e[0] / np.sqrt(1 - phi**2)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x + mean


def ar1_fast(n: int, phi: float, seed: int) -> np.ndarray:
    from scipy.signal import lfilter

    e = np.random.default_rng(seed).standard_normal(n)
    return np.asarray(lfilter([1.0], [1.0, -phi], e))


def test_autocorrelation_matches_direct() -> None:
    x = ar1(300, 0.6, 1)
    rho = autocorrelation(x)
    d = x - x.mean()
    direct = np.array([np.sum(d[: len(d) - t] * d[t:]) for t in range(len(d))]) / np.sum(d * d)
    np.testing.assert_allclose(rho, direct, atol=1e-10)
    assert rho[0] == 1.0
    c = autocorrelation(np.full(10, 3.0))
    assert c[0] == 1.0 and np.all(c[1:] == 0.0)


def test_tau_int_ar1_phi09() -> None:
    phi = 0.9
    exact = 0.5 * (1 + phi) / (1 - phi)  # 9.5
    x = ar1_fast(200_000, phi, 7)
    tau = integrated_autocorr_time(x)
    assert tau == pytest.approx(exact, rel=0.10)
    assert effective_sample_size(x) == pytest.approx(len(x) / (2 * exact), rel=0.10)


def test_tau_int_iid_is_half() -> None:
    x = np.random.default_rng(0).standard_normal(50_000)
    assert integrated_autocorr_time(x) == pytest.approx(0.5, abs=0.08)
    assert integrated_autocorr_time(np.ones(100)) == 0.5


def test_split_rhat_iid_and_shifted() -> None:
    rng = np.random.default_rng(1)
    iid = [rng.standard_normal(5000) for _ in range(4)]
    assert split_rhat(iid) <= 1.01
    shifted = [rng.standard_normal(5000) + 0.5 * k for k in range(4)]
    assert split_rhat(shifted) > 1.1
    assert split_rhat([np.ones(10), np.ones(10)]) == 1.0
    # una sola cadena con deriva: el split lo detecta
    drift = np.linspace(0, 5, 4000) + rng.standard_normal(4000)
    assert split_rhat([drift, drift.copy()]) > 1.1
    with pytest.raises(ValueError):
        split_rhat([np.ones(3)])


def test_geweke_stationary_vs_drift() -> None:
    x = ar1_fast(20_000, 0.5, 3)
    assert abs(geweke_z(x)) < 3.0
    drift = x + np.linspace(0, 3, x.size)
    assert abs(geweke_z(drift)) > 3.0
    with pytest.raises(ValueError):
        geweke_z(x[:5])


def test_post_burn_in_pure() -> None:
    x = np.arange(10.0)
    y = post_burn_in(x, 0.5)
    np.testing.assert_array_equal(y, np.arange(5.0, 10.0))
    y[0] = -1.0
    assert x[5] == 5.0
    with pytest.raises(ValueError):
        post_burn_in(x, 1.0)


def _chain(samples: dict[str, np.ndarray]) -> ChainResult:
    n = len(next(iter(samples.values())))
    return ChainResult(
        engine=Engine.METROPOLIS, theta=0.1, samples=samples,
        sample_steps=np.arange(1, n + 1, dtype=np.int64), w_final=np.zeros((2, 2)),
        states=(), acceptance=0.4, step_size=0.1, n_steps=n,
    )


def test_ensemble_summary_equilibrated_and_not() -> None:
    cfg = EnsembleConfig()
    good = [
        _chain({"action": ar1_fast(4000, 0.8, s) + 10.0, "mean_weight": 0.5 + 0.01 * ar1_fast(4000, 0.8, s + 50)})
        for s in range(4)
    ]
    s = ensemble_summary(good, ["action", "mean_weight"], cfg, 0.5)
    assert s.equilibrated and s.geweke_ok
    assert s.means["action"] == pytest.approx(10.0, abs=0.5)
    exact_tau = 0.5 * 1.8 / 0.2
    assert s.tau_int["action"] == pytest.approx(exact_tau, rel=0.35)
    assert s.ess["action"] == pytest.approx(4 * 2000 / (2 * exact_tau), rel=0.35)
    assert s.rhat["action"] <= 1.05 and s.ses["action"] > 0
    # cadenas con medias distintas -> no equilibrado
    bad = [_chain({"action": ar1_fast(4000, 0.8, s) + 5.0 * s, "mean_weight": np.full(4000, 0.5) + 0.0 * s})
           for s in range(4)]
    sb = ensemble_summary(bad, ["action", "mean_weight"], cfg, 0.5)
    assert not sb.equilibrated and sb.rhat["action"] > 1.05
    # ESS insuficiente
    short = [_chain({"action": ar1_fast(60, 0.9, s)}) for s in range(4)]
    ss = ensemble_summary(short, ["action"], cfg, 0.5)
    assert not ss.equilibrated
    with pytest.raises(ValueError):
        ensemble_summary(good[:1], ["action"], cfg, 0.5)
