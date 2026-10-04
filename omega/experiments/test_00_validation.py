"""Experimento 0 (M§22, ANALYSIS §6.3): validacion del algoritmo de geometria sobre redes conocidas.

Si falla una fila obligatoria no se ejecutan los Experimentos 1-5 (M§22). Las tolerancias son las
del analisis; no se ensanchan. Semillas explicitas via SeedSequence (D-28).
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
import pytest

from omega.config.settings import DimensionConfig, GraphConfig, SpectralConfig
from omega.experiments.reference_graphs import (
    complete_graph,
    open_lattice,
    periodic_lattice,
    random_geometric_torus,
    random_regular,
)
from omega.geometry.observables import geometry_observables
from omega.network.initialization import random_uniform_weights
from omega.network.weights import permute
from omega.types import FloatArray, GeometryObservables

DIM = DimensionConfig()
SPEC = SpectralConfig()
ENTROPY = 20240601
EXP_ID = 0


def _rng(*key: int) -> np.random.Generator:
    """Generator(PCG64(SeedSequence(ENTROPY, spawn_key=(EXP_ID, *key))))."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=(EXP_ID, *key))))


def _geo(w: FloatArray, w_min: float = 0.1) -> GeometryObservables:
    return geometry_observables(w, GraphConfig(w_min=w_min), DIM, SPEC)


# (nombre, forma, D_eff +- tol, D_s +- tol) segun §6.3
PERIODIC: tuple[tuple[str, tuple[int, ...], float, float, float, float], ...] = (
    ("anillo300", (300,), 1.0, 0.05, 1.0, 0.10),
    ("anillo1000", (1000,), 1.0, 0.05, 1.0, 0.10),
    ("toro17x17", (17, 17), 2.0, 0.10, 2.0, 0.15),
    ("toro32x32", (32, 32), 2.0, 0.10, 2.0, 0.15),
    ("toro7^3", (7, 7, 7), 3.0, 0.30, 3.0, 0.25),
    ("toro10^3", (10, 10, 10), 3.0, 0.25, 3.0, 0.25),
    ("toro15^3", (15, 15, 15), 3.0, 0.15, 3.0, 0.25),
)


@lru_cache(maxsize=None)
def _periodic_geo(shape: tuple[int, ...]) -> GeometryObservables:
    return _geo(periodic_lattice(shape))


@pytest.mark.parametrize(("name", "shape", "d_eff", "tol_d", "d_s", "tol_s"), PERIODIC, ids=[p[0] for p in PERIODIC])
def test_periodic_lattices(name: str, shape: tuple[int, ...], d_eff: float, tol_d: float, d_s: float, tol_s: float) -> None:
    g = _periodic_geo(shape)
    assert g.d_eff.status == "ok", name
    assert g.d_eff.plateau, name
    assert abs(g.d_eff.value - d_eff) <= tol_d, f"{name}: D_eff={g.d_eff.value:.4f} (esperado {d_eff}+-{tol_d})"
    assert g.d_eff_hops.value == pytest.approx(g.d_eff.value, abs=1e-6)  # W=1: distancia = saltos
    assert g.d_s.status == "ok", name
    assert abs(g.d_s.value - d_s) <= tol_s, f"{name}: D_s={g.d_s.value:.4f} (esperado {d_s}+-{tol_s})"


# (dim, grado medio, D_eff esperado, tol); metrica euclidiana, N=1000, 4 semillas
RGG: tuple[tuple[int, float, float, float], ...] = ((2, 20.0, 2.0, 0.15), (3, 40.0, 3.0, 0.30))


@pytest.mark.parametrize(("dim", "degree", "expected", "tol"), RGG, ids=["rgg2d", "rgg3d"])
def test_random_geometric_torus(dim: int, degree: float, expected: float, tol: float) -> None:
    values = []
    for rep in range(4):
        w = random_geometric_torus(1000, dim, degree, _rng(1, dim, rep))
        g = _geo(w)
        assert g.d_eff.status == "ok"
        values.append(g.d_eff.value)
    assert abs(float(np.mean(values)) - expected) <= tol, f"RGG {dim}D: D_eff={np.mean(values):.4f} {values}"
    assert all(abs(v - expected) <= tol for v in values), f"RGG {dim}D por semilla: {values}"


@pytest.mark.parametrize("n", [100, 300])
@pytest.mark.parametrize("w_min", [0.0, 0.5, 0.9])
@pytest.mark.parametrize("rep", [0, 1])
def test_negative_control_uniform_weights(n: int, w_min: float, rep: int) -> None:
    w = random_uniform_weights(n, _rng(2, n, int(w_min * 10), rep))
    g = _geo(w, w_min)
    assert g.d_eff.status == "no_window", f"U(0,1) N={n} W_min={w_min}: D_eff={g.d_eff.value}"
    assert math.isnan(g.d_eff.value)


@pytest.mark.parametrize("n", [100, 300])
def test_negative_control_complete_graph(n: int) -> None:
    g = _geo(complete_graph(n))
    assert g.d_eff.status == "no_window" and g.d_eff_ball.status == "no_window"
    assert g.d_eff_hops.status == "no_window" and g.d_s.status == "no_window"


def test_negative_control_random_3_regular() -> None:
    """3-regular: no_window o plateau=False (ANALYSIS §6.3)."""
    for rep in range(4):
        g = _geo(random_regular(300, 3, _rng(3, rep)))
        assert g.d_eff.status == "no_window" or not g.d_eff.plateau, f"3-regular rep={rep}: {g.d_eff}"


def _same(a: float, b: float, tol: float = 1e-9) -> bool:
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tol


def _fields(g: GeometryObservables) -> list[tuple[str, float]]:
    out = [("L", g.path_length), ("L_hops", g.path_length_hops), ("diam", g.diameter)]
    for tag, est in (("d_eff", g.d_eff), ("ball", g.d_eff_ball), ("hops", g.d_eff_hops), ("d_s", g.d_s)):
        out += [(tag, est.value), (tag + "_se", est.stderr)]
    return out


def _permutation_cases() -> list[tuple[str, FloatArray, float]]:
    return [
        ("U(0,1) N=100", random_uniform_weights(100, _rng(4, 0)), 0.9),
        ("U(0,1) N=100 Wmin=0.5", random_uniform_weights(100, _rng(4, 1)), 0.5),
        ("toro 12x12", periodic_lattice((12, 12)), 0.1),
        ("RGG 2D", random_geometric_torus(300, 2, 14.0, _rng(4, 2)), 0.1),
        ("3-regular", random_regular(120, 3, _rng(4, 3)), 0.1),
    ]


@pytest.mark.parametrize("case", range(5))
def test_permutation_invariance(case: int) -> None:
    name, w, w_min = _permutation_cases()[case]
    perm = _rng(5, case).permutation(w.shape[0]).astype(np.int64)
    a, b = _geo(w, w_min), _geo(permute(w, perm), w_min)
    for (tag, x), (_, y) in zip(_fields(a), _fields(b), strict=True):
        assert _same(x, y), f"{name}: {tag} {x} vs {y}"
    for ea, eb in ((a.d_eff, b.d_eff), (a.d_eff_ball, b.d_eff_ball), (a.d_eff_hops, b.d_eff_hops), (a.d_s, b.d_s)):
        assert ea.status == eb.status and ea.window == eb.window


@pytest.mark.parametrize("shape", [(300,), (17, 17), (7, 7, 7)], ids=["anillo", "toro2d", "toro3d"])
def test_wmin_independence_for_unit_weights(shape: tuple[int, ...]) -> None:
    ref = _periodic_geo(shape)
    w = periodic_lattice(shape)
    for w_min in (0.0, 0.01, 0.5, 0.9, 0.99):
        g = _geo(w, w_min)
        for (tag, x), (_, y) in zip(_fields(ref), _fields(g), strict=True):
            assert _same(x, y), f"W_min={w_min}: {tag} {x} vs {y}"
        assert g.d_eff.status == ref.d_eff.status and g.d_s.status == ref.d_s.status


def test_open_lattice_is_report_only() -> None:
    """Redes abiertas y estimador bola: solo informe (no hay tolerancia obligatoria), pero deben estimarse."""
    g = _geo(open_lattice((7, 7, 7)))
    assert g.d_eff.status == "ok" and g.d_eff_ball.status == "ok"
    assert g.d_eff_ball.value < g.d_eff.value  # la bola esta mas sesgada (R3)
