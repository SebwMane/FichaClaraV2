"""Tests del modulo omega.causality (bloque C0, Rev. 3)."""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from omega.causality.estimators import (
    between_counts,
    bd_action_2d,
    interval_abundances,
    midpoint_dimension,
    mm_dimension,
    mm_fraction,
    ordering_fraction,
)
from omega.causality.exact import expected_abundance_2order, expected_bd_action_2order
from omega.causality.orders import (
    is_partial_order,
    kr_order,
    random_2order,
    transitive_closure,
    transitive_percolation,
)
from omega.causality.sprinkling import causal_matrix, diamond_volume, interval_volume, sprinkle_diamond


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _chain(n: int) -> np.ndarray:
    return np.triu(np.ones((n, n), dtype=bool), k=1)


def test_mm_fraction_known_values() -> None:
    assert mm_fraction(1) == pytest.approx(1.0)
    assert mm_fraction(2) == pytest.approx(0.5)
    assert mm_fraction(4) == pytest.approx(0.1)


def test_interval_volume_2d_and_4d() -> None:
    assert interval_volume(1.0, 2) == pytest.approx(0.5)
    assert interval_volume(1.0, 4) == pytest.approx(math.pi / 24)
    assert diamond_volume(2) == pytest.approx(0.5)


def test_chain_and_antichain() -> None:
    assert mm_dimension(_chain(10)) == pytest.approx(1.0)
    assert math.isnan(mm_dimension(np.zeros((10, 10), dtype=bool)))


def test_bd_action_small_cases() -> None:
    # cadena de 3: N=3, N0=2 links, N1=1 -> S = 3 - 4 + 4 = 3
    assert bd_action_2d(_chain(3)) == 3
    # antichain de 5: S = 5
    assert bd_action_2d(np.zeros((5, 5), dtype=bool)) == 5
    # cadena de 4: N0=3, N1=2, N2=1 -> 4 - 6 + 8 - 2 = 4
    assert bd_action_2d(_chain(4)) == 4


def test_between_counts_marks_unrelated() -> None:
    b = between_counts(_chain(4))
    assert b[0, 3] == 2 and b[3, 0] == -1 and b[0, 1] == 0
    assert list(interval_abundances(_chain(4), 3)) == [3, 2, 1, 0]


def test_sprinkling_inside_diamond_and_partial_order() -> None:
    for d in (2, 3, 4):
        x = sprinkle_diamond(200, d, _rng(d))
        r = np.linalg.norm(x[:, 1:], axis=1)
        assert np.all(r <= 0.5 - np.abs(x[:, 0]) + 1e-12)
        assert is_partial_order(causal_matrix(x))


def test_random_2order_matches_lightcone_sprinkling() -> None:
    c = random_2order(300, _rng(1))
    assert is_partial_order(c)
    assert ordering_fraction(c) == pytest.approx(0.5, abs=0.05)


def test_controls_are_partial_orders_and_layered() -> None:
    k = kr_order(200, _rng(2))
    assert is_partial_order(k)
    p = transitive_percolation(200, 0.05, _rng(3))
    assert is_partial_order(p)
    # AUDITORIA C0 (A-3): el test original era vacuo ('or True') y su afirmacion era falsa
    # (un par capa1-capa3 tiene ~n/8 elementos entre ellos). Propiedad correcta: altura <= 3,
    # es decir, no existe cadena de 4 elementos (C^3 = 0).
    kf = k.astype(np.float32)
    assert not np.any(kf @ kf @ kf)
    assert int(between_counts(k).max()) > 0


def test_transitive_closure_simple() -> None:
    adj = np.zeros((3, 3), dtype=bool)
    adj[0, 1] = adj[1, 2] = True
    assert transitive_closure(adj)[0, 2]


def test_reproducible() -> None:
    a = sprinkle_diamond(50, 3, _rng(9))
    b = sprinkle_diamond(50, 3, _rng(9))
    assert np.array_equal(a, b)


def test_midpoint_dimension_2d_reasonable() -> None:
    c = random_2order(512, _rng(4))
    assert 1.5 < midpoint_dimension(c) < 2.6


def test_exact_bd_mean_is_one() -> None:
    for n in (5, 40, 300):
        assert expected_bd_action_2order(n) == pytest.approx(1.0, abs=1e-6)


def test_abundances_match_exact_expectation() -> None:
    n, reps = 128, 200
    ab = np.array([interval_abundances(random_2order(n, _rng(500 + i)), 2) for i in range(reps)])
    for k in range(3):
        se = ab[:, k].std(ddof=1) / math.sqrt(reps)
        assert abs(ab[:, k].mean() - expected_abundance_2order(k, n)) <= 4 * se


def test_thick_layer_is_partial_order() -> None:
    from omega.causality.orders import thick_layer_order

    c = thick_layer_order(160, 8, 0.5, _rng(11))
    assert is_partial_order(c)


def test_profile_helpers() -> None:
    from omega.causality.profile import chain_abundances, d_ref, profile_distance, reference_abundances

    assert list(chain_abundances(5, 3)) == [4, 3, 2, 1]
    # una cadena coincide exactamente con su referencia
    assert profile_distance(interval_abundances(_chain(30), 8), reference_abundances(1, 30)) == 0.0
    assert d_ref(2.4) == 2 and d_ref(0.3) == 1 and d_ref(7.0) == 4 and d_ref(float("nan")) is None
    ref2 = reference_abundances(2, 64)
    assert ref2[0] == pytest.approx(expected_abundance_2order(0, 64))


def test_layered_2order_action_formula() -> None:
    from omega.causality.orders import layered_2order

    n, k = 64, 4
    c = layered_2order(n, k, _rng(12))
    assert is_partial_order(c)
    m = n // k
    # links solo entre capas consecutivas; capas no consecutivas tienen >= m >= 3 elementos entre ellas
    assert bd_action_2d(c) == n - 2 * (k - 1) * m * m


def test_bipartite_all_relations_are_links() -> None:
    from omega.causality.orders import random_bipartite

    c = random_bipartite(20, 0.5, _rng(13))
    assert is_partial_order(c)
    assert bd_action_2d(c) == 20 - 2 * int(c.sum())


def test_phase_helpers() -> None:
    from omega.causality.phase import bipartite_phase_exact, phase_average

    mag, se = phase_average([2, 2, 2, 2], 0.7)
    assert mag == pytest.approx(1.0) and se == pytest.approx(0.0, abs=1e-12)
    assert bipartite_phase_exact(8, 0.0) == 1.0
    assert bipartite_phase_exact(8, math.pi / 2) == pytest.approx(0.0, abs=1e-12)


def test_smeared_action_limit_and_exact_mean() -> None:
    from omega.causality.estimators import bd_action_2d_smeared
    from omega.causality.exact import expected_smeared_action_2order

    c = random_2order(40, _rng(21))
    assert bd_action_2d_smeared(c, 1 - 1e-9) == pytest.approx(bd_action_2d(c), abs=1e-3)
    assert expected_smeared_action_2order(30, 1 - 1e-9) == pytest.approx(1.0, abs=1e-3)


def test_mcmc_moves_keep_permutation_and_beta0_accepts_all() -> None:
    from omega.causality.mcmc import run_chain

    v0 = _rng(22).permutation(12).astype(np.int64)
    res = run_chain(v0, 0.0, 2, 4, 2, _rng(23))
    assert res.acceptance == 1.0
    assert len(res.samples) == 2
    for c in res.samples:
        assert is_partial_order(c)


def test_wanglandau_chunking_is_deterministic_and_n4_exact() -> None:
    from omega.causality.wanglandau import BATCH, exact_dos, new_state, wl_run

    d = exact_dos(4)
    assert sum(d.values()) == 24
    a = new_state(4, 4 - 8, 4 + 8, _rng(31))
    b = new_state(4, 4 - 8, 4 + 8, _rng(31))
    a = wl_run(a, 2 * BATCH)
    b = wl_run(wl_run(b, BATCH), BATCH)
    assert a.lng == b.lng and a.v == b.v and a.t == b.t == 2 * BATCH


def test_bd_action_4d_small_cases() -> None:
    from omega.causality.estimators import bd_action_4d, bd_action_4d_smeared

    # cadena de 5: N0=4, N1=3, N2=2, N3=1 -> 5 - 4 + 27 - 32 + 8 = 4
    assert bd_action_4d(_chain(5)) == 4
    assert bd_action_4d(np.zeros((6, 6), dtype=bool)) == 6
    c = causal_matrix(sprinkle_diamond(60, 4, _rng(41)))
    assert bd_action_4d_smeared(c, 1 - 1e-9) == pytest.approx(bd_action_4d(c), abs=1e-3)


def test_slab_torus_relation_and_retarded_sum() -> None:
    from omega.causality.estimators import bd_action_4d, bd_action_4d_smeared
    from omega.causality.torus import retarded_action_4d, slab_causal_matrix, sprinkle_slab

    # dos puntos separados espacialmente 0.9 en x estan a 0.1 por la imagen minima
    x = np.array([[0.0, 0.05, 0.5, 0.5], [0.15, 0.95, 0.5, 0.5]])
    assert slab_causal_matrix(x)[0, 1]
    c = slab_causal_matrix(sprinkle_slab(150, _rng(51)))
    assert is_partial_order(c)
    assert retarded_action_4d(c, 1.0).sum() == pytest.approx(bd_action_4d(c))
    assert retarded_action_4d(c, 0.25).sum() == pytest.approx(bd_action_4d_smeared(c, 0.25))


def test_general_mcmc_consistency() -> None:
    from omega.causality.estimators import bd_action_4d, between_counts
    from omega.causality.mcmc_general import init_state, sweep

    x = sprinkle_diamond(24, 4, _rng(61))
    x = x[np.argsort(x[:, 0])]
    st = init_state(causal_matrix(x))
    assert st.s == pytest.approx(bd_action_4d(st.c))
    r = _rng(62)
    for _ in range(20):
        sweep(st, 0.05, r)
    assert is_partial_order(st.c)
    full = between_counts(st.c)
    assert np.array_equal(st.b[st.c], full[st.c])
    assert st.s == pytest.approx(bd_action_4d(st.c))


def test_ssee_full_region_is_zero_and_pairs() -> None:
    from omega.causality.entropy import pauli_jordan, sj_split, ssee

    c = random_2order(120, _rng(71))
    idelta = pauli_jordan(c)
    idt, w = sj_split(idelta, 1e-9)
    allidx = np.arange(120)
    assert abs(ssee(idt, w, allidx, truncate=False)) < 1e-6  # vacio puro en la region completa
    sub = allidx[:40]
    s = ssee(idt, w, sub, truncate=False)
    assert s > 0


def test_interval_closest_rule() -> None:
    from omega.causality.entropy import interval_closest

    r = interval_closest(_chain(10), 4)
    assert r.size == 4 and list(r) == [0, 1, 2, 3]


def test_e1prime_helpers() -> None:
    from omega.causality.entropy import causal_complement, interval_balanced, knee_cutoff, pauli_jordan

    c = _chain(12)
    r = interval_balanced(c, 4)
    assert 3 <= r.size <= 5  # la ventana N/4 +- 5% se redondea hacia afuera (floor/ceil)
    assert causal_complement(c, r).size == 0  # en una cadena todo esta relacionado
    anti = np.zeros((6, 6), dtype=bool)
    assert interval_balanced(anti, 3).size == 0
    k = knee_cutoff(pauli_jordan(random_2order(200, _rng(81))))
    assert k is None or k > 0


def test_knee_powerlaw_synthetic() -> None:
    from omega.causality.entropy import knee_powerlaw

    # espectro sintetico: 1/k hasta el rango 40 y luego cae mucho mas rapido
    k = np.arange(1, 201, dtype=float)
    lam = np.where(k <= 40, 100 / k, 100 / 40 * (40 / k) ** 10)
    m = np.diag(np.concatenate([lam, -lam])).astype(np.complex128)
    thr, beta = knee_powerlaw(m, 1600)
    assert beta == pytest.approx(-1.0, abs=1e-6)
    # rango 41: cae 20% (bajo el umbral de 30%); rango 42: cae 36% -> rodilla en lambda_42 = 2.5 (40/42)^10
    assert thr == pytest.approx(2.5 * (40 / 42) ** 10, rel=1e-9)


def test_sj_from_eig_and_knee_from_spectrum_match() -> None:
    from omega.causality.entropy import knee_from_spectrum, knee_powerlaw, pauli_jordan, sj_both, sj_from_eig

    c = random_2order(150, _rng(91))
    idelta = pauli_jordan(c)
    mu, u = np.linalg.eigh(idelta)
    _, idt_a, w_a = sj_both(idelta, 0.7)
    idt_b, w_b = sj_from_eig(mu, u, 0.7)
    assert np.allclose(idt_a, idt_b) and np.allclose(w_a, w_b)
    assert knee_from_spectrum(mu, 150)[1] == pytest.approx(knee_powerlaw(idelta, 150)[1])


def _f5prima() -> ModuleType:
    """tools/f5prima.py (F5′); se importa por ruta igual que la ejecuta el orquestador."""
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("f5prima")


def test_f5prima_polyhedral_cones() -> None:
    f = _f5prima()
    t_half = 0.5
    for norm, seed in (("max", 401), ("l1", 402)):
        x = f.sprinkle_polyhedral(300, norm, _rng(seed))
        assert x.shape == (300, 4)
        nx = np.max(np.abs(x[:, 1:]), axis=1) if norm == "max" else np.sum(np.abs(x[:, 1:]), axis=1)
        assert np.all(np.abs(x[:, 0]) <= t_half) and np.all(nx <= t_half - np.abs(x[:, 0]))
        c = f.polyhedral_matrix(x, norm)
        assert is_partial_order(c) and np.count_nonzero(c) > 0
        assert np.array_equal(c, f.polyhedral_cone(norm)(300, _rng(seed)))
    # definicion: octaedrica ⊆ Minkowski ⊆ cubica sobre los mismos puntos
    x = sprinkle_diamond(300, 4, _rng(403))
    mink, cub, octa = causal_matrix(x), f.polyhedral_matrix(x, "max"), f.polyhedral_matrix(x, "l1")
    assert not np.any(mink & ~cub) and not np.any(octa & ~mink)
    assert np.count_nonzero(cub) > np.count_nonzero(mink) > np.count_nonzero(octa) > 0
    with pytest.raises(ValueError):
        f.polyhedral_cone("l2")


def test_f5prima_de_sitter_and_frw() -> None:
    f = _f5prima()
    de_sitter_4d = importlib.import_module("f3f4").de_sitter_4d
    assert np.array_equal(f.de_sitter(4)(250, _rng(411)), de_sitter_4d(250, _rng(411)))
    for d in (2, 3, 4, 5):
        for pts, gen in ((f.de_sitter_points, f.de_sitter(d)), (f.frw_points, f.frw(d))):
            x = pts(250, d, _rng(420 + d))
            assert x.shape == (250, d)
            r = np.linalg.norm(x[:, 1:], axis=1)
            assert np.all(r <= 0.5 - np.abs(x[:, 0]))
            c = gen(250, _rng(420 + d))
            assert np.array_equal(c, causal_matrix(x)) and is_partial_order(c)
    # ambas densidades crecen con t: de Sitter (0.5/|t - 1|)^d y FRW ((t + 1)/1.5)^(2d)
    assert np.mean(f.de_sitter_points(400, 4, _rng(431))[:, 0]) > 0.05
    assert np.mean(f.frw_points(400, 4, _rng(432))[:, 0]) > 0.05


def _f7() -> ModuleType:
    """tools/f7.py (F7); se importa por ruta igual que _f5prima."""
    tools = str(Path(__file__).resolve().parents[1] / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("f7")


def test_f7_lp_cone() -> None:
    f, g = _f7(), _f5prima()
    t_half = 0.5
    for p, seed in ((1.5, 501), (3.0, 502), (6.0, 503)):
        x = f.sprinkle_lp(300, p, _rng(seed))
        assert x.shape == (300, 4)
        nx = np.sum(np.abs(x[:, 1:]) ** p, axis=1) ** (1 / p)
        assert np.all(np.abs(x[:, 0]) <= t_half) and np.all(nx <= t_half - np.abs(x[:, 0]))
        c = f.lp_matrix(x, p)
        assert is_partial_order(c) and np.count_nonzero(c) > 0
        assert np.array_equal(c, f.lp_cone(p)(300, _rng(seed)))
    # sobre los mismos puntos: octaedrica ⊆ L^1.5 ⊆ L^3 ⊆ L^6 ⊆ cubica, con mas relaciones al crecer p
    x = sprinkle_diamond(300, 4, _rng(504))
    chain = [g.polyhedral_matrix(x, "l1")] + [f.lp_matrix(x, p) for p in (1.5, 3.0, 6.0)] + [g.polyhedral_matrix(x, "max")]
    for lo, hi in zip(chain, chain[1:]):
        assert not np.any(lo & ~hi)
        assert np.count_nonzero(hi) > np.count_nonzero(lo)
    for bad in (0.5, math.inf):
        with pytest.raises(ValueError):
            f.lp_cone(bad)


def test_f7_de_sitter_strong_and_radiation() -> None:
    f, g = _f7(), _f5prima()
    for d in (3, 4):
        # la version estandar reproduce exactamente a f5prima.de_sitter(d) y frw_power_points(2d) a frw_points
        assert np.array_equal(f.de_sitter_eta(d)(250, _rng(511 + d)), g.de_sitter(d)(250, _rng(511 + d)))
        assert np.array_equal(f.frw_power_points(250, d, 2.0 * d, _rng(521 + d)), g.frw_points(250, d, _rng(521 + d)))
    cases = (
        (f.de_sitter_eta_points(250, 4, _rng(531), -1.2, -0.2), f.de_sitter_eta(4, -1.2, -0.2)(250, _rng(531))),
        (f.frw_power_points(250, 4, 4.0, _rng(532)), f.frw_radiation(4)(250, _rng(532))),
    )
    for x, c in cases:
        assert x.shape == (250, 4)
        assert np.all(np.linalg.norm(x[:, 1:], axis=1) <= 0.5 - np.abs(x[:, 0]))
        assert np.array_equal(c, causal_matrix(x)) and is_partial_order(c)
        assert np.mean(x[:, 0]) > 0.02  # ambas densidades crecen con t (medias ≈ 0.12 y 0.06)
    # mas curvatura: el de Sitter fuerte se concentra mas arriba que el estandar
    strong, std = f.de_sitter_eta_points(600, 4, _rng(533), -1.2, -0.2), f.de_sitter_eta_points(600, 4, _rng(533))
    assert np.mean(strong[:, 0]) > np.mean(std[:, 0])
    with pytest.raises(ValueError):
        f.de_sitter_eta_points(10, 4, _rng(534), -1.5, 0.5)


def test_f7_shape_distance_shift_invariant() -> None:
    f = _f7()
    ref = np.array([100.0, 80, 60, 45, 30, 20, 12, 7, 3])
    nk = np.array([120.0, 70, 66, 40, 33, 18, 15, 5, 4])
    v = f.log_ratio(nk, ref)
    assert f.shape_distance(nk, ref) == pytest.approx(float(np.max(np.abs(v - v.mean()))))
    for c in (-0.7, 0.3, 2.0):
        # sumar c a todos los v_k equivale a N_k + 1 -> (N_k + 1) e^c
        shifted = (nk + 1) * math.exp(c) - 1
        assert np.allclose(f.log_ratio(shifted, ref), v + c)
        assert f.shape_distance(shifted, ref) == pytest.approx(f.shape_distance(nk, ref), abs=1e-12)
    assert f.shape_distance((ref + 1) * 1.7 - 1, ref) == pytest.approx(0.0, abs=1e-12)
    assert len(f.TEST) == 17 and list(f.TEST).index("cono_p6_4d") == 10 and f.TEST["cono_p6_4d"][0] == "informativa"
