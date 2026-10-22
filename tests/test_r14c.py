"""Pruebas rápidas de tools/r14c.py (R14-C v2). Semillas propias (880001), nunca la entropía del experimento."""

from __future__ import annotations

import importlib
import itertools
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
from scipy.stats import energy_distance

from omega.causality.estimators import bd_action_4d, between_counts, mm_dimension
from omega.causality.mcmc_general import height
from omega.causality.orders import transitive_percolation
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond


def _m() -> ModuleType:
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("r14c")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(880001, spawn_key=(seed,)))


def _brute_chains(c: np.ndarray, k: int) -> int:
    n = c.shape[0]
    tot = 0
    for sub in itertools.permutations(range(n), k):
        if all(c[sub[i], sub[i + 1]] for i in range(k - 1)):
            tot += 1
    return tot


def test_chain_counts_brute_force() -> None:
    m = _m()
    orders = [causal_matrix(sprinkle_diamond(9, 4, _rng(1))), causal_matrix(sprinkle_diamond(10, 2, _rng(2))),
              transitive_percolation(9, 0.3, _rng(3))]
    perm = _rng(4).permutation(9)
    orders.append(orders[2][np.ix_(perm, perm)])  # etiquetado no natural
    for c in orders:
        counts, _ = m.chain_counts(c)
        for k in range(2, 6):
            exp = _brute_chains(c, k)
            got = counts[k - 1] if k <= len(counts) else 0
            assert got == exp, (k, got, exp)
        assert counts[0] == c.shape[0]


def test_longest_chain_matches_height() -> None:
    m = _m()
    for s in range(6):
        c = transitive_percolation(60, 0.08, _rng(10 + s))
        nat = np.argsort(c.sum(axis=0), kind="stable")  # número de predecesores: extensión lineal
        assert m.longest_chain(c) == height(c[np.ix_(nat, nat)])
        x = sprinkle_diamond(80, 4, _rng(20 + s))
        x = x[np.argsort(x[:, 0])]  # orden temporal = etiquetado natural
        c2 = causal_matrix(x)
        assert m.longest_chain(c2) == height(c2)
        perm = _rng(30 + s).permutation(80)
        assert m.longest_chain(c2[np.ix_(perm, perm)]) == height(c2)
    assert m.longest_chain(np.zeros((5, 5), dtype=bool)) == 1
    assert m.longest_chain(np.triu(np.ones((7, 7), dtype=bool), 1)) == 7
    assert m.chain_counts(np.zeros((4, 4), dtype=bool), kmax=5)[0] == [4, 0, 0, 0, 0]


def test_expected_chains_vs_monte_carlo() -> None:
    """⟨C_k⟩ = n(n−1)···(n−k+1)·χ_k con n fijo (factorial descendente), frente a Monte Carlo en diamantes planos 4D."""
    m = _m()
    n, reps = 40, 1500
    ck = np.array([m.chain_counts(causal_matrix(sprinkle_diamond(n, 4, _rng(1000 + s))), kmax=5)[0][1:5] for s in range(reps)],
                  dtype=np.float64)
    for i, k in enumerate((2, 3, 4, 5)):
        th = m.expected_chains(n, k)
        z = (ck[:, i].mean() - th) / (ck[:, i].std(ddof=1) / math.sqrt(reps))
        assert abs(z) < 4, (k, z)
    assert math.isclose(m.expected_chains(n, 2), n * (n - 1) / 20)
    # χ_k = 12^(k−1)·6 / [k (2k−1)! (2k+1)!]
    for k in (2, 3, 4, 5):
        chi = 12 ** (k - 1) * 6 / (k * math.factorial(2 * k - 1) * math.factorial(2 * k + 1))
        assert math.isclose(m.CHI[k], chi, rel_tol=1e-12)


def test_interval_measures_conventions() -> None:
    m = _m()
    mm = 30
    x = sprinkle_diamond(mm, 4, _rng(40))
    tips = np.array([[-0.5, 0, 0, 0], [0.5, 0, 0, 0]])
    ct = causal_matrix(np.vstack([tips[:1], x, tips[1:]]))
    inner = np.arange(1, mm + 1)
    cin = ct[np.ix_(inner, inner)]
    b = between_counts(ct)
    res = m.interval_measures(cin, b[np.ix_(inner, inner)], ct)
    assert math.isclose(res["r"], np.count_nonzero(cin) / (mm * (mm - 1) / 2))
    assert int(res["N_k"].sum()) == res["R"] == np.count_nonzero(cin)
    assert np.array_equal(res["N_k"][:5], np.bincount(between_counts(cin)[cin], minlength=5)[:5])  # B global = B del interior
    assert math.isclose(res["s_media"], bd_action_4d(cin) / mm)
    assert res["C_k"][0] == res["R"]
    d_r = mm_dimension(cin)
    assert math.isfinite(res["dmid_menos_dr"]) and math.isfinite(d_r)


def _order_from_links(n: int, edges: list[tuple[int, int]]) -> np.ndarray:
    c = np.zeros((n, n), dtype=bool)
    for i, j in edges:
        c[i, j] = True
    for k in range(n):  # cierre transitivo
        c |= c[:, [k]] & c[[k], :]
    return c


def _links(c: np.ndarray) -> np.ndarray:
    out: np.ndarray = between_counts(c) == 0
    return out


def test_laplacian_path_bipartite_isolated() -> None:
    m = _m()
    # cadena de 6 elementos: Hasse = camino P6; espectro del Laplaciano normalizado 1 − cos(πk/5)
    c = np.triu(np.ones((6, 6), dtype=bool), 1)
    w = m.laplacian_spectrum(_links(c))
    assert np.allclose(np.sort(w), np.sort(1 - np.cos(np.pi * np.arange(6) / 5)), atol=1e-10)
    assert math.isclose(w[-1], 2.0, abs_tol=1e-10)  # el camino es bipartito
    # orden de dos capas completo (3 abajo, 4 arriba): K_{3,4}, espectro {0, 1 (×5), 2}
    c = np.zeros((7, 7), dtype=bool)
    c[:3, 3:] = True
    w = m.laplacian_spectrum(_links(c))
    assert np.allclose(w, [0, 1, 1, 1, 1, 1, 2], atol=1e-10)
    s = m.spectrum_summary(w)
    assert math.isclose(s["lambda2"], 1.0) and math.isclose(s["lambda_max"], 2.0) and sum(s["hist"]) == 7
    # vértice aislado: 𝒟^(−1/2) = 0 (pseudoinversa) -> autovalor 1; el resto no cambia
    c2 = np.zeros((8, 8), dtype=bool)
    c2[:7, :7] = c
    lap = m.normalized_laplacian(_links(c2))
    assert lap[7, 7] == 1.0 and not np.any(lap[7, :7]) and not np.any(lap[:7, 7])
    w2 = m.laplacian_spectrum(_links(c2))
    assert np.allclose(w2, [0, 1, 1, 1, 1, 1, 1, 2], atol=1e-10)
    # grados: d_in + d_out de los links
    d = m.link_degrees(_links(_order_from_links(4, [(0, 1), (1, 2), (0, 3)])))
    assert d.tolist() == [2, 2, 1, 1]


def test_laplacian_spectrum_in_range_on_sprinkling() -> None:
    m = _m()
    c = causal_matrix(sprinkle_diamond(300, 4, _rng(50)))
    w = m.laplacian_spectrum(_links(c))
    assert w.size == 300 and w[0] > -1e-9 and w[-1] < 2 + 1e-9
    assert abs(w[0]) < 1e-9  # componente principal


def test_energy_distance_properties() -> None:
    m = _m()
    r = _rng(60)
    a, b = r.normal(size=200), r.normal(size=150) + 0.5
    assert energy_distance(a, a) == 0.0
    assert math.isclose(energy_distance(a, b), energy_distance(b, a))
    assert energy_distance(a, b) > energy_distance(a, r.normal(size=150))
    # forma V: D² = 2E|X − Y| − E|X − X'| − E|Y − Y'|
    exy = np.abs(a[:, None] - b[None, :]).mean()
    exx = np.abs(a[:, None] - a[None, :]).mean()
    eyy = np.abs(b[:, None] - b[None, :]).mean()
    assert math.isclose(energy_distance(a, b) ** 2, 2 * exy - exx - eyy, rel_tol=1e-9)
    # con pesos = muestras repetidas (distribuciones discretas: grados, abundancias, histogramas)
    xs = np.arange(5, dtype=np.float64)
    w1, w2 = np.array([1, 3, 0, 2, 4.0]), np.array([2, 2, 2, 2, 2.0])
    rep1, rep2 = np.repeat(xs, w1.astype(int)), np.repeat(xs, w2.astype(int))
    assert math.isclose(m.ed_weighted(xs, w1, xs, w2), energy_distance(rep1, rep2), rel_tol=1e-12)
    assert m.ed_weighted(xs, w1, xs, 3 * w1) == 0.0  # solo cuenta la distribución, no el total
    assert math.isnan(m.ed_weighted(xs, np.zeros(5), xs, w2))


def test_holm() -> None:
    m = _m()
    p = [0.01, 0.04, 0.03, 0.005, math.nan]
    adj = m.holm(p)
    # orden 0.005, 0.01, 0.03, 0.04 con k = 4: 0.02, 0.03, max(0.03, 0.06), max(0.06, 0.04)
    assert np.allclose(adj[:4], [0.03, 0.06, 0.06, 0.02])
    assert math.isnan(adj[4])
    assert m.holm([0.9, 0.8]) == [1.0, 1.0]
    q = list(_rng(70).random(30))
    a = m.holm(q)
    order = np.argsort(q)
    assert all(a[order[i]] <= a[order[i + 1]] for i in range(29))
    assert all(a[i] >= q[i] for i in range(30))


def test_family_effect_and_slope() -> None:
    m = _m()
    r = _rng(80)
    null = m.family_effect([r.normal(size=10) for _ in range(3)], _rng(81))
    shifted = m.family_effect([r.normal(size=10) + 3 for _ in range(3)], _rng(81))
    assert shifted["p"] < 1e-6 < null["p"]
    assert shifted["ic95"][0] < shifted["delta"] < shifted["ic95"][1]
    assert m.family_effect([], _rng(82))["n_estructuras"] == 0
    row = {"por_m": {str(mm): {"medible": True, "r_media": 0.1 + 0.2 * math.sqrt(mm / 4000)} for mm in m.MS}}
    assert math.isclose(m.slope_a(row, 4000), 0.2, rel_tol=1e-9)
    row["por_m"]["40"]["medible"] = row["por_m"]["80"]["medible"] = row["por_m"]["160"]["medible"] = False
    assert math.isnan(m.slope_a(row, 4000))


def test_measure_small_and_mc() -> None:
    m = _m()
    c = causal_matrix(sprinkle_diamond(900, 4, _rng(90)))
    row = m.measure(c, _rng(91))
    assert row["N"] == 900 and sum(row["grado_bincount"]) == 900
    assert row["links"] == sum(i * v for i, v in enumerate(row["grado_bincount"])) // 2
    cell = row["por_m"]["20"]
    assert cell["medible"] and cell["n_sel"] == len(cell["r"]) <= 300
    assert 0.07 < cell["r_media"] < 0.13
    assert math.isclose(sum(cell["abund_media"]), 1.0)  # media de N_k/R sobre intervalos con R > 0
    assert 0 <= cell["abund_sin_relaciones"] < 10
    assert not row["por_m"]["160"]["medible"]
    assert 0 < row["laplaciano"]["lambda2"] < row["laplaciano"]["lambda_max"] <= 2 + 1e-9
    mc = m.mc_rows(20, 30)
    assert mc["reps"] == 30 and len(mc["L"]) == 30 and len(mc["r"]) == 30
    # determinismo de la referencia de Monte Carlo
    assert m.mc_rows(20, 5)["r"] == mc["r"][:5]
