"""WP-D: observables de transicion sobre datos sinteticos."""

from __future__ import annotations

import numpy as np
import pytest

from omega.statistics.transition import (
    bimodality_coefficient,
    binder_cumulant,
    fluctuation,
    order_parameter_histogram,
    susceptibility,
    transition_summary,

)


def _rng(seed: int) -> np.random.Generator:
    """Generador PCG64 explícito (prohibido np.random global)."""
    return np.random.Generator(np.random.PCG64(seed))


def test_fluctuation_and_susceptibility() -> None:
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert fluctuation(x) == pytest.approx(1.25)
    assert susceptibility(x, 10) == pytest.approx(12.5)
    with pytest.raises(ValueError):
        susceptibility(x, 0)


def test_binder_known_distributions() -> None:
    rng = _rng(0)
    assert binder_cumulant(rng.standard_normal(400_000)) == pytest.approx(0.0, abs=0.02)
    assert binder_cumulant(rng.uniform(size=400_000)) == pytest.approx(0.4, abs=0.02)
    assert binder_cumulant(np.tile([-1.0, 1.0], 500)) == pytest.approx(2.0 / 3.0, abs=1e-12)
    assert np.isnan(binder_cumulant(np.ones(10)))


def test_bimodality_coefficient() -> None:
    rng = _rng(1)
    uni = rng.normal(0.5, 0.05, 50_000)
    bi = np.concatenate([rng.normal(0.2, 0.03, 25_000), rng.normal(0.8, 0.03, 25_000)])
    assert bimodality_coefficient(bi) > 5.0 / 9.0 > bimodality_coefficient(uni)
    assert bimodality_coefficient(uni) == pytest.approx(1.0 / 3.0, abs=0.03)
    assert bimodality_coefficient(rng.uniform(size=400_000)) == pytest.approx(5.0 / 9.0, abs=0.01)
    assert np.isnan(bimodality_coefficient(np.zeros(5)))


def test_histogram_density_integrates_to_one() -> None:
    m = _rng(2).uniform(size=5000)
    dens, edges = order_parameter_histogram(m, 20)
    assert dens.shape == (20,) and edges.shape == (21,)
    assert float(np.sum(dens * np.diff(edges))) == pytest.approx(1.0)


def test_transition_summary_bimodal_peak() -> None:
    rng = _rng(3)
    lams = [0.0, 0.25, 0.5, 0.75, 1.0]
    samples = {}
    for lam in lams:
        if lam == 0.5:  # coexistencia de fases
            m = np.concatenate([rng.normal(0.2, 0.01, 5000), rng.normal(0.8, 0.01, 5000)])
        else:
            m = rng.normal(0.2 + 0.6 * lam, 0.01, 10_000)
        samples[lam] = m
    # insercion desordenada: debe ordenar
    shuffled = {k: samples[k] for k in [0.75, 0.0, 1.0, 0.5, 0.25]}
    ts = transition_summary(shuffled, n_edges=28)
    np.testing.assert_array_equal(ts.lambdas, lams)
    assert ts.peak_lambda == 0.5
    assert ts.susceptibility[2] == ts.susceptibility.max()
    assert ts.bimodality[2] > 5.0 / 9.0
    assert all(ts.bimodality[i] < 5.0 / 9.0 for i in (0, 1, 3, 4))
    assert ts.binder[2] == pytest.approx(2.0 / 3.0, abs=0.02)
    assert ts.mean[1] == pytest.approx(0.35, abs=0.001)
    assert ts.derivative[0] == pytest.approx((ts.mean[1] - ts.mean[0]) / 0.25)
    with pytest.raises(ValueError):
        transition_summary({0.0: samples[0.0]}, 10)


def test_transition_derivative_linear_central() -> None:
    lams = np.linspace(0, 1, 6)
    samples = {float(l): np.array([2.0 * l + 1.0, 2.0 * l + 1.0, 2.0 * l + 1.0]) for l in lams}
    ts = transition_summary(samples, 3)
    np.testing.assert_allclose(ts.derivative, 2.0, atol=1e-12)
