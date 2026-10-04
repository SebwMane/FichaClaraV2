"""Pruebas de omega.geometry.weyl (DESIGN §1.3, §3.2 WP-B)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from omega.config.settings11 import WeylConfig
from omega.experiments.reference_graphs import (
    complete_graph,
    periodic_lattice,
    random_geometric_torus,
)
from omega.geometry.weyl import fiedler_length, laplacian_spectrum, weyl_dimension, weyl_staircase
from omega.types import FloatArray

CFG = WeylConfig()


@pytest.mark.parametrize(
    ("shape", "expected", "tol"),
    [
        ((300,), 1.01, 0.05),
        ((17, 17), 2.13, 0.2),
        ((28, 28), 2.08, 0.2),
        ((7, 7, 7), 2.88, 0.3),
        ((9, 9, 9), 2.98, 0.25),
    ],
)
def test_weyl_table_fast(shape: tuple[int, ...], expected: float, tol: float) -> None:
    e = weyl_dimension(periodic_lattice(shape), 0.5, CFG)
    assert e.status == "ok" and e.plateau
    assert e.value == pytest.approx(expected, abs=tol)
    assert e.method == "weyl_combinatorial_binary"


@pytest.mark.slow
@pytest.mark.parametrize(
    ("shape", "expected", "tol"),
    [((1000,), 1.01, 0.05), ((32, 32), 2.07, 0.2), ((10, 10, 10), 3.08, 0.25), ((15, 15, 15), 3.15, 0.25)],
)
def test_weyl_table_slow(shape: tuple[int, ...], expected: float, tol: float) -> None:
    e = weyl_dimension(periodic_lattice(shape), 0.5, CFG)
    assert e.status == "ok" and e.value == pytest.approx(expected, abs=tol)


@pytest.mark.slow
def test_weyl_rgg() -> None:
    d3 = [
        weyl_dimension(random_geometric_torus(800, 3, 12, np.random.Generator(np.random.PCG64(s)), euclidean=False), 0.5, CFG).value
        for s in range(3)
    ]
    assert float(np.mean(d3)) == pytest.approx(2.96, abs=0.25)
    for s in range(3):
        w = random_geometric_torus(800, 2, 10, np.random.Generator(np.random.PCG64(s)), euclidean=False)
        assert weyl_dimension(w, 0.5, CFG).value == pytest.approx(2.0, abs=0.25)


def test_complete_graph_no_window() -> None:
    e = weyl_dimension(complete_graph(40), 0.5, CFG)
    assert e.status == "no_window" and math.isnan(e.value)


def test_tiny_component_insufficient() -> None:
    e = weyl_dimension(periodic_lattice((8,)), 0.5, CFG)
    assert e.status == "insufficient_component"


def _uniform(n: int, seed: int) -> FloatArray:
    u = np.triu(np.random.Generator(np.random.PCG64(seed)).random((n, n)), 1)
    return np.asarray(u + u.T, dtype=np.float64)


@pytest.mark.parametrize("w_min", [0.1, 0.5, 0.9])
def test_uniform_nulls_dimension_above_five(w_min: float) -> None:
    e = weyl_dimension(_uniform(100, 1), w_min, CFG)
    assert e.status == "ok" and e.value > 5.0


def test_strict_threshold_edges() -> None:
    w = periodic_lattice((100,)) * 0.5
    assert weyl_dimension(w, 0.5, CFG).status == "insufficient_component"  # W == w_min: sin aristas
    assert weyl_dimension(w, 0.49, CFG).status == "ok"


def test_staircase_monotone_without_zero_mode() -> None:
    eigs = laplacian_spectrum(periodic_lattice((30,)), "combinatorial")
    assert abs(eigs[0]) < 1e-9
    levels, counts = weyl_staircase(eigs, 1e-9)
    assert np.all(levels > 1e-9) and np.all(np.diff(levels) > 0) and np.all(np.diff(counts) > 0)
    assert counts[-1] == 29.0  # n-1 autovalores no nulos
    assert counts[0] == 2.0  # degeneracion 2 del anillo agrupada en un nivel
    lv, ct = weyl_staircase(laplacian_spectrum(complete_graph(6), "combinatorial"), 1e-9)
    assert lv.size == 1 and ct[0] == 5.0


def test_normalized_spectrum_range() -> None:
    e = laplacian_spectrum(periodic_lattice((6, 6)), "normalized")
    assert abs(e[0]) < 1e-9 and e[-1] <= 2.0 + 1e-9
    cfg = WeylConfig(laplacian="normalized")
    assert weyl_dimension(periodic_lattice((17, 17)), 0.5, cfg).method == "weyl_normalized_binary"


def test_fiedler_length() -> None:
    n = 100
    lam2 = 2.0 - 2.0 * math.cos(2.0 * math.pi / n)
    assert fiedler_length(periodic_lattice((n,)), 0.5) == pytest.approx(lam2**-0.5)
    assert math.isnan(fiedler_length(np.zeros((4, 4)), 0.5))


def test_b18_weyl_and_fiedler_equal_separate_calls() -> None:
    from omega.geometry.weyl import weyl_and_fiedler

    rng = np.random.Generator(np.random.PCG64(3))
    graphs = [periodic_lattice((17, 17)), periodic_lattice((60,)), np.zeros((4, 4)), periodic_lattice((6, 6))]
    u = np.triu(rng.random((80, 80)), 1)
    graphs.append(np.asarray(u + u.T, dtype=np.float64))
    for cfg in (WeylConfig(), WeylConfig(laplacian="normalized"), WeylConfig(graph="thresholded_weighted")):
        for w in graphs:
            est, xi = weyl_and_fiedler(w, 0.5, cfg)
            ref = weyl_dimension(w, 0.5, cfg)
            ref_xi = fiedler_length(w, 0.5)
            assert est.method == ref.method and est.status == ref.status and est.plateau == ref.plateau and est.window == ref.window
            assert np.array_equal(est.value, ref.value, equal_nan=True) and np.array_equal(est.stderr, ref.stderr, equal_nan=True)
            assert np.array_equal(est.scales, ref.scales) and np.array_equal(est.profile, ref.profile)
            assert np.array_equal(est.local_slopes, ref.local_slopes, equal_nan=True)
            assert (math.isnan(xi) and math.isnan(ref_xi)) or xi == ref_xi


def test_b18_single_spectrum_computation(monkeypatch: pytest.MonkeyPatch) -> None:
    import omega.geometry.weyl as mod

    calls: list[str] = []
    real = mod.laplacian_spectrum

    def spy(w: FloatArray, kind: str) -> FloatArray:
        calls.append(kind)
        return real(w, kind)  # type: ignore[arg-type]

    monkeypatch.setattr(mod, "laplacian_spectrum", spy)
    mod.weyl_and_fiedler(periodic_lattice((17, 17)), 0.5, WeylConfig())
    assert calls == ["combinatorial"]
