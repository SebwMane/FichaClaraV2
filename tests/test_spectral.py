"""Pruebas de omega.geometry.spectral y observables (ANALYSIS §6.2: spectral)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.settings import DimensionConfig, GraphConfig, SpectralConfig
from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.geometry.observables import geometry_observables
from omega.geometry.spectral import (
    heat_spectrum,
    mean_return_probability,
    spectral_dimension,
    time_grid,
    walk_operator_spectrum,
)
from omega.network.initialization import random_uniform_weights

SP = SpectralConfig()


@pytest.mark.parametrize("q", [0.25, 0.5, 1.0])
def test_lazy_spectrum_in_range(q: float) -> None:
    mu = walk_operator_spectrum(periodic_lattice((8, 8)), q)
    assert mu.min() >= 1.0 - 2.0 * q - 1e-12 and mu.max() <= 1.0 + 1e-12
    assert mu.max() == pytest.approx(1.0, abs=1e-12)


def test_lazy_spectrum_random_weights_in_range() -> None:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(21)))
    mu = walk_operator_spectrum(random_uniform_weights(50, rng), 0.5)
    assert mu.min() >= -1e-12 and mu.max() <= 1.0 + 1e-12


def test_return_probability_decreasing_to_one_over_n() -> None:
    w = periodic_lattice((6, 6))
    mu = walk_operator_spectrum(w, 0.5)
    t = np.array([0.0, 1.0, 2.0, 5.0, 20.0, 100.0, 1e4])
    p = mean_return_probability(np.maximum(mu, 0.0), t, "lazy_walk")
    assert p[0] == pytest.approx(1.0)
    assert np.all(np.diff(p) <= 1e-15)
    assert p[-1] == pytest.approx(1.0 / 36.0, rel=1e-9)
    ph = mean_return_probability(heat_spectrum(w), np.array([0.0, 1.0, 10.0, 1e4]), "heat_normalized")
    assert ph[0] == pytest.approx(1.0) and np.all(np.diff(ph) <= 1e-15)
    assert ph[-1] == pytest.approx(1.0 / 36.0, rel=1e-9)


def test_return_probability_matches_matrix_power() -> None:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(22)))
    w = random_uniform_weights(15, rng)
    p_mat = 0.5 * np.eye(15) + 0.5 * w / w.sum(axis=1, keepdims=True)
    direct = np.trace(np.linalg.matrix_power(p_mat, 7)) / 15
    got = mean_return_probability(walk_operator_spectrum(w, 0.5), np.array([7.0]), "lazy_walk")[0]
    assert got == pytest.approx(direct, rel=1e-10)


def test_time_grid() -> None:
    t = time_grid(100, SP)
    assert t[0] == 1.0 and t[-1] == 1000.0 and np.all(np.diff(t) > 0) and np.all(t == np.rint(t))


def test_complete_graph_no_window() -> None:
    est = spectral_dimension(complete_graph(100), SP)
    assert est.status == "no_window" and np.isnan(est.value) and est.window is None


def test_ring_and_torus_values() -> None:
    ring = spectral_dimension(periodic_lattice((300,)), SP)
    assert ring.status == "ok" and ring.value == pytest.approx(1.0, abs=0.1)
    t2 = spectral_dimension(periodic_lattice((17, 17)), SP)
    assert t2.status == "ok" and t2.value == pytest.approx(2.0, abs=0.15)
    assert t2.window is not None
    lo, hi = t2.window
    assert t2.scales[hi] / t2.scales[lo] >= 3.0 and t2.scales[lo] >= SP.t_min
    assert np.all(t2.profile[lo : hi + 1] >= SP.saturation_factor / 289)


def test_heat_method_runs_and_insufficient() -> None:
    est = spectral_dimension(periodic_lattice((200,)), SpectralConfig(method="heat_normalized"))
    assert est.status in ("ok", "no_window")
    assert spectral_dimension(np.zeros((20, 20)), SP).status == "insufficient_component"


def test_geometry_observables_ring() -> None:
    o = geometry_observables(periodic_lattice((300,)), GraphConfig(), DimensionConfig(), SP)
    assert o.d_eff.value == pytest.approx(1.0, abs=1e-9) and o.d_eff_hops.value == pytest.approx(1.0, abs=1e-9)
    assert o.d_s.status == "ok"
    assert o.diameter == pytest.approx(150.0, rel=1e-6)
    assert o.path_length_hops == pytest.approx(300**2 / 4 / 299, rel=1e-9)


def test_geometry_observables_empty_graph() -> None:
    w = np.zeros((20, 20))
    o = geometry_observables(w, GraphConfig(), DimensionConfig(), SP)
    assert o.d_eff.status == "insufficient_component" and o.d_s.status == "insufficient_component"
    assert np.isnan(o.path_length) and np.isnan(o.diameter)
