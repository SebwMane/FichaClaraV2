"""Tests del modulo omega.causality (bloque C0, Rev. 3)."""

from __future__ import annotations

import math

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
    assert r.size == 4
    assert causal_complement(c, r).size == 0  # en una cadena todo esta relacionado
    anti = np.zeros((6, 6), dtype=bool)
    assert interval_balanced(anti, 3).size == 0
    k = knee_cutoff(pauli_jordan(random_2order(200, _rng(81))))
    assert k is None or k > 0
