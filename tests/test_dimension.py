"""Pruebas de omega.geometry.dimension (ANALYSIS §6.2: dimension)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.settings import DimensionConfig
from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.geometry.dimension import (
    effective_dimension,
    fit_loglog,
    local_slopes,
    profile_dimension,
    radius_grid,
    scaling_window,
    shell_profile,
)
from omega.geometry.distances import distance_matrix, hop_distance_matrix, threshold_adjacency

BALL = DimensionConfig(estimator="ball")
SHELL = DimensionConfig(estimator="shell")
HUGE = 10**12


@pytest.mark.parametrize("dim", [1.0, 2.5, 3.0])
def test_power_law_ball_recovers_dimension(dim: float) -> None:
    r = np.geomspace(1.0, 1000.0, 60)
    est = profile_dimension(r, 2.0 * r**dim, HUGE, BALL, integer=False)
    assert est.status == "ok" and est.plateau
    assert est.value == pytest.approx(dim, abs=1e-6)


@pytest.mark.parametrize("dim", [1.0, 2.5, 3.0])
def test_power_law_shell_recovers_dimension(dim: float) -> None:
    r = np.geomspace(1.0, 1000.0, 6000)
    est = profile_dimension(r, 2.0 * r**dim, HUGE, SHELL, integer=False)
    assert est.status == "ok"
    assert est.value == pytest.approx(dim, abs=1e-6)


@pytest.mark.parametrize("est", [BALL, SHELL])
def test_step_profile_gives_no_window(est: DimensionConfig) -> None:
    r = np.arange(1.0, 21.0)
    prof = np.where(r < 5.0, 1.0, 1000.0)
    out = profile_dimension(r, prof, 1000, est, integer=True)
    assert out.status == "no_window" and np.isnan(out.value) and out.window is None


def test_scaling_window_rules() -> None:
    s = np.arange(1.0, 11.0)
    p = np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1000.0])
    assert scaling_window(s, p, 100.0, 3.0, 3) == (0, 5)
    assert scaling_window(s, p, 10.0, 3.0, 3) == (0, 2)  # razon exactamente 3, 3 puntos
    assert scaling_window(s, p, 10.0, 3.0, 4) is None  # solo 3 puntos <= 10
    assert scaling_window(s, p, 9.0, 4.0, 3) is None
    assert scaling_window(s, p, 100.0, 7.0, 3) is None  # razon 6 < 7


def test_fit_and_local_slopes() -> None:
    x = np.geomspace(1.0, 100.0, 20)
    slope, err, r2 = fit_loglog(x, 3.0 * x**1.7)
    assert slope == pytest.approx(1.7, abs=1e-12) and err < 1e-12 and r2 == pytest.approx(1.0)
    assert np.allclose(local_slopes(x, 3.0 * x**1.7), 1.7)
    with pytest.raises(ValueError):
        fit_loglog(x, -x)


def test_shell_profile_integer_and_continuous() -> None:
    s, v = shell_profile(np.array([1.0, 2.0, 3.0]), np.array([5.0, 13.0, 25.0]), True)
    assert s.tolist() == [1.0, 2.0, 3.0] and v.tolist() == [4.0, 8.0, 12.0]
    s2, v2 = shell_profile(np.array([1.0, 4.0]), np.array([2.0, 11.0]), False)
    assert s2.tolist() == [2.0] and v2.tolist() == [3.0]


def test_radius_grid_integer_and_continuous() -> None:
    d = hop_distance_matrix(threshold_adjacency(periodic_lattice((20,)), 0.5))
    assert radius_grid(d, SHELL).tolist() == [float(i) for i in range(1, 11)]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(5)))
    w = periodic_lattice((20,)) * rng.uniform(0.3, 1.0, (20, 20))
    w = np.triu(w, 1)
    w = w + w.T
    dc = distance_matrix(w, 0.0, 1e-9)
    g = radius_grid(dc, SHELL)
    assert g.size == SHELL.n_radii and np.all(np.diff(g) > 0) and g[-1] == pytest.approx(dc.max())


@pytest.mark.parametrize(("shape", "expected", "tol"), [((300,), 1.0, 1e-9), ((20, 20), 2.0, 1e-9)])
def test_lattice_shell_exact_low_dim(shape: tuple[int, ...], expected: float, tol: float) -> None:
    d = distance_matrix(periodic_lattice(shape), 0.5, 1e-9)
    est = effective_dimension(d, SHELL)
    assert est.status == "ok" and est.plateau
    assert est.value == pytest.approx(expected, abs=tol)


def test_complete_graph_no_window() -> None:
    d = distance_matrix(complete_graph(80), 0.0, 1e-9)
    est = effective_dimension(d, SHELL)
    assert est.status == "no_window" and np.isnan(est.value)
    assert effective_dimension(d, BALL).status == "no_window"


def test_insufficient_component_and_disconnected() -> None:
    w = np.zeros((30, 30))
    w[0, 1] = w[1, 0] = 1.0
    d = distance_matrix(w, 0.5, 1e-9)
    assert effective_dimension(d, SHELL).status == "insufficient_component"
    # dos anillos: se usa solo la gigante (mayor)
    big = periodic_lattice((60,))
    small = periodic_lattice((20,))
    w2 = np.zeros((80, 80))
    w2[:60, :60] = big
    w2[60:, 60:] = small
    est = effective_dimension(distance_matrix(w2, 0.5, 1e-9), SHELL)
    assert est.status == "ok" and est.value == pytest.approx(1.0, abs=1e-9)


def test_trim_first_literal_differs_on_3d() -> None:
    d = distance_matrix(periodic_lattice((10, 10, 10)), 0.5, 1e-9)
    assert effective_dimension(d, SHELL, integer_trim=1).value > effective_dimension(d, SHELL, integer_trim=0).value
