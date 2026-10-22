"""Pruebas rápidas de tools/r14b.py (R14-B v2). Semillas propias (870001), nunca la entropía del experimento."""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np

from omega.causality.estimators import interval_abundances, ordering_fraction
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond


def _m() -> ModuleType:
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("r14b")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(870001, spawn_key=(seed,)))


def test_dilution_is_induced_order() -> None:
    m = _m()
    x = sprinkle_diamond(400, 4, _rng(1))
    c = causal_matrix(x)
    for q in (0.5, 0.25, 0.125):
        idx, cd = m.dilute(c, q, _rng(2))
        assert np.all(np.diff(idx) > 0)
        assert np.array_equal(cd, causal_matrix(x[idx]))  # orden inducido = orden de los puntos conservados
        # mismo generador -> misma selección; la selección es un uniforme por elemento
        u = _rng(2).random(400)
        assert np.array_equal(idx, np.flatnonzero(u < q))
    idx1, c1 = m.dilute(c, 1.0, None)
    assert idx1.size == 400 and c1 is c
    # tamaño ~ Binomial(n, q)
    sizes = [m.dilute(c, 0.25, _rng(10 + s))[0].size for s in range(200)]
    assert abs(np.mean(sizes) - 100) < 3


def test_twins_on_constructed_order() -> None:
    m = _m()
    # 0 < 1 < 2, 0 < 3 < 2 (3 gemelo de 1), 4 y 5 aislados (también gemelos entre sí), 6 sobre 2.
    n = 7
    c = np.zeros((n, n), dtype=np.bool_)
    for i, j in [(0, 1), (1, 2), (0, 3), (3, 2), (0, 2), (2, 6), (0, 6), (1, 6), (3, 6)]:
        c[i, j] = True
    row = m.measure(c, _rng(3))
    assert row["pares_gemelos"] == 2
    assert math.isclose(row["frac_gemelos"], 2 / 21)
    c[:, 6] = False  # 6 queda aislado: gemelo de 4 y 5 (3 pares); 1 y 3 siguen siendo gemelos (1 par)
    assert m.measure(c, _rng(3))["pares_gemelos"] == 1 + 3


def test_measure_matches_estimators() -> None:
    m = _m()
    c = causal_matrix(sprinkle_diamond(600, 4, _rng(4)))
    row = m.measure(c, _rng(5))
    pairs = 600 * 599 / 2
    assert math.isclose(row["frac_rel"], ordering_fraction(c))
    assert row["N_k"] == [int(v) for v in interval_abundances(c, 8)]
    assert np.allclose(row["perfil"], np.array(row["N_k"]) / pairs)
    for mm in m.M_R:
        cell = row["por_m"][str(mm)]
        assert cell["local"] == (mm <= 0.1 * row["m_max"])
        if cell["medible"]:
            assert 30 <= cell["n_sel"] <= 300
    v = m.obs_vector(row)
    assert v.shape == (len(m.COMPONENTS),)


def test_ledoit_wolf_and_band_coverage() -> None:
    """Datos gaussianos sintéticos (p = 14, casi colineales). El p95 de los 40 de ajuste subcubre (~74% con p = 14); la
    banda final (19.º de 20 de validación) cubre en promedio 19/21 ≈ 0.905 de puntos nuevos (§3)."""
    m = _m()
    r = _rng(6)
    p = 14
    a = r.normal(size=(p, p))
    cov = a @ a.T + 0.05 * np.eye(p)  # casi colineal
    chol = np.linalg.cholesky(cov)
    floors = np.full(p, 1e-12)
    m_comps = tuple(range(p))
    covs = []
    for t in range(60):
        data = r.normal(size=(60, p)) @ chol.T
        band = m.Band(data[:40], data[40:], m_comps, floors)
        new = r.normal(size=(500, p)) @ chol.T
        covs.append(float(np.mean(band.dist(new, sub=True) <= band.radio)))
        assert band.radio >= band.r_fit and 0 <= band.lam <= 1
    cov_mean = float(np.mean(covs))
    assert 0.88 <= cov_mean <= 0.99, cov_mean
    s, lam = m.ledoit_wolf(r.normal(size=(40, 5)))
    assert np.all(np.linalg.eigvalsh(s) > 0) and 0 <= lam <= 1


def test_twin_floor_band() -> None:
    """Una componente constante en la referencia (√gemelos = 0) no hace singular la covarianza y un par de gemelos
    cuenta como unas 2 sd (piso 0.5)."""
    m = _m()
    r = _rng(7)
    data = np.column_stack([r.normal(size=60), np.zeros(60)])
    floors = np.array([1e-12, 0.5])
    band = m.Band(data[:40], data[40:], (0, 1), floors)
    d0 = band.dist(np.array([[0.0, 0.0]]), sub=True)[0]
    d1 = band.dist(np.array([[0.0, 1.0]]), sub=True)[0]
    assert math.isfinite(d1) and 1.5 < d1 - d0 < 2.5


def test_labels() -> None:
    m = _m()
    assert m.state(0.2, 0.9) == "dentro" and m.state(1.1, 3.0) == "fuera" and m.state(0.8, 1.2) == "ambiguo"
    lab = m.label
    assert lab(["dentro", "fuera", "fuera", "fuera"], False) == "en banda microscópica"
    assert lab(["fuera", "fuera", "dentro", "dentro"], False) == "converge a banda"
    assert lab(["fuera", "dentro", "dentro", "dentro"], True) == "no converge"  # crece de 1/2 a 1/8
    assert lab(["ambiguo", "fuera", "dentro", "dentro"], False) == "converge a banda"
    assert lab(["fuera", "fuera", "fuera", "dentro"], False) == "no converge"
    assert lab(["fuera", "fuera", "ambiguo", "dentro"], False) == "indeterminada"
    assert lab(["ambiguo", "ambiguo", "ambiguo", "ambiguo"], False) == "indeterminada"
    # trayectoria completa: D/radio que baja de 4 a 0.3 con poco ruido -> converge
    r = _rng(8)
    dn = np.array([4.0, 2.0, 0.5, 0.3])[None, :] * np.exp(0.05 * r.normal(size=(10, 4)))
    tr = m.trajectory(dn, _rng(9), n_boot=300)
    assert tr["etiqueta"] == "converge a banda" and tr["estados"] == ["fuera", "fuera", "dentro", "dentro"]
    assert not tr["cambio_1_8_menos_1_2"]["crece"]
    tr2 = m.trajectory(0.4 * np.exp(0.05 * r.normal(size=(10, 4))), _rng(10), n_boot=300)
    assert tr2["etiqueta"] == "en banda microscópica" and tr2["en_banda_todo_q"]
    # dentro en 1/4 y 1/8 pero creciendo claramente desde 1/2 -> no converge
    tr3 = m.trajectory(np.array([3.0, 0.3, 0.5, 0.8])[None, :] * np.exp(0.03 * r.normal(size=(10, 4))), _rng(11), n_boot=300)
    assert tr3["cambio_1_8_menos_1_2"]["crece"] and tr3["etiqueta"] == "no converge"


def test_structure_table() -> None:
    m = _m()
    assert len(m.STRUCTS) == 17 and len(m.GEOMS) == 8
    assert m.primary_target("D5_gemelos") == "plana" and m.primary_target("C0_kasner") == "C0_kasner"
    # cada generador produce un orden con n elementos (n chico)
    for name in m.NAMES:
        c, _ = m.GEN[name](120, _rng(11), _rng(12))
        assert c.shape == (120, 120) and not np.any(c & c.T)
