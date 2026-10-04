"""Evidencia por corrida (omega.certificate.evidence): controles, igualdad con Omega-1.0, determinismo y tiempo."""

from __future__ import annotations

import dataclasses
import json
import time

import numpy as np

from omega.certificate.evidence import (
    _cached_topology_summary,
    collect_run_evidence,
    evidence_rng,
    evidence_summary,
    weight_stats,
)
from omega.certificate.taxonomy import assess_run
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import Omega11Config
from omega.contracts import FailureCode
from omega.controls.random_geometric import rgg_torus
from omega.controls.small_world import watts_strogatz
from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.phases.scan import default_config, observe
from omega.topology.persistence import topology_summary
from omega.types import DimensionEstimate, FloatArray, RunStatus

KEY = SeedKey(20240901, (7, 0, 0))


def _cfg(n: int) -> Omega11Config:
    return Omega11Config(base=default_config(n, master_entropy=20240901, replicates=3))


def _same_estimate(a: DimensionEstimate, b: DimensionEstimate) -> bool:
    return (
        a.status == b.status
        and a.window == b.window
        and a.plateau == b.plateau
        and a.method == b.method
        and (a.value == b.value or (np.isnan(a.value) and np.isnan(b.value)))
        and np.array_equal(a.scales, b.scales)
        and np.array_equal(a.profile, b.profile, equal_nan=True)
    )


def test_rgg3_binary_n800_passes_and_time() -> None:
    cfg = _cfg(800)
    w = rgg_torus(800, 3, 12.0, make_rng(KEY))
    t0 = time.perf_counter()
    ev = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY))
    elapsed = time.perf_counter() - t0
    a = assess_run(ev, cfg)
    assert a.passes, (a.codes, [k for k, v in a.flags.items() if not v])
    assert ev.dimension_class == 3
    assert elapsed < 20.0, f"evidencia N=800 tardo {elapsed:.1f} s"


def test_watts_strogatz_p005_primary_f3() -> None:
    cfg = _cfg(800)
    w = watts_strogatz(800, 12, 0.05, make_rng(KEY))
    a = assess_run(collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY)), cfg)
    assert a.primary is FailureCode.F3 and not a.passes
    assert a.flags["locality"] is False and a.flags["d_volume"] is False


def test_geometry_equals_omega10_observe() -> None:
    cfg = _cfg(150)
    w = rgg_torus(150, 2, 10.0, make_rng(KEY))
    ev = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY))
    obs = observe(w, cfg.base)
    assert ev.topology == obs.topology
    g0, g1 = ev.geometry, obs.geometry
    for name in ("d_eff", "d_eff_ball", "d_eff_hops", "d_s"):
        assert _same_estimate(getattr(g0, name), getattr(g1, name)), name
    for name in ("path_length", "path_length_hops", "diameter"):
        assert getattr(g0, name) == getattr(g1, name), name


def test_determinism_by_seed_key() -> None:
    cfg = _cfg(150)
    w = rgg_torus(150, 3, 12.0, make_rng(KEY))
    a = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY), with_qrc=True)
    b = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY), with_qrc=True)
    assert evidence_summary(a) == evidence_summary(b)
    assert a.isotropy == b.isotropy and a.annulus == b.annulus
    assert np.array_equal(a.curvature.edge_values, b.curvature.edge_values)
    assert a.qrc is not None and b.qrc is not None
    assert np.array_equal(a.qrc.ratio_mean, b.qrc.ratio_mean, equal_nan=True)


def test_evidence_rng_is_independent_stream() -> None:
    r1 = evidence_rng(KEY).random(4)
    assert np.array_equal(r1, make_rng(SeedKey(KEY.entropy, KEY.spawn_key + (2000,))).random(4))
    assert not np.array_equal(r1, make_rng(KEY).random(4))
    assert isinstance(evidence_rng(KEY).bit_generator, np.random.PCG64)


def test_qrc_only_on_request() -> None:
    cfg = _cfg(100)
    w = periodic_lattice((10, 10))
    assert collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY)).qrc is None
    ev = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY), with_qrc=True)
    assert ev.qrc is not None and ev.qrc.deltas.shape[0] == len(cfg.curvature.qrc_deltas)


def test_weight_stats() -> None:
    w = np.zeros((4, 4))
    assert weight_stats(w) == (0.0, 0.0)
    assert weight_stats(complete_graph(5)) == (1.0, 0.0)
    m, cv = weight_stats(periodic_lattice((6,)))
    assert abs(m - 6 / 15) < 1e-12 and cv > 0.0


def test_trivial_states_are_robust() -> None:
    cfg = _cfg(30)
    empty = collect_run_evidence(np.zeros((30, 30)), RunStatus.CONVERGED, cfg, evidence_rng(KEY))
    assert assess_run(empty, cfg).primary is FailureCode.F0
    full = collect_run_evidence(complete_graph(30), RunStatus.CONVERGED, cfg, evidence_rng(KEY))
    assert assess_run(full, cfg).primary is FailureCode.F1
    nonconv = collect_run_evidence(complete_graph(30), RunStatus.MAX_STEPS, cfg, evidence_rng(KEY))
    assert assess_run(nonconv, cfg).primary is FailureCode.F10
    tiny = collect_run_evidence(complete_graph(3), RunStatus.CONVERGED, _cfg(3), evidence_rng(KEY))
    assert tiny.n == 3


def test_summary_is_strict_json() -> None:
    cfg = _cfg(60)
    for w in (np.zeros((60, 60)), periodic_lattice((60,)), complete_graph(60)):
        s = evidence_summary(collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(KEY)))
        json.dumps(s, allow_nan=False)


def test_cached_topology_equals_wp_c_summary() -> None:
    cfg = _cfg(120)
    rng = make_rng(KEY)
    states: list[FloatArray] = [rgg_torus(120, 2, 8.0, rng)]
    u = np.triu(rng.random((120, 120)), 1)
    u = np.where(u > 0.9, u, 0.0)
    states.append(np.asarray(u + u.T, dtype=np.float64))  # no binario: niveles distintos por theta
    for w in states:
        args = (w, 0.1, cfg.base.graph.w_min_sensitivity, cfg.certificate.g_connected, cfg.topology)
        a, b = topology_summary(*args), _cached_topology_summary(*args)
        assert a.stable_flags == b.stable_flags and a.stable == b.stable
        assert a.at_w_min == b.at_w_min and a.clique_at_w_min == b.clique_at_w_min
        for f in dataclasses.fields(a.curves):
            va, vb = getattr(a.curves, f.name), getattr(b.curves, f.name)
            assert (va == vb) if isinstance(va, tuple) else np.array_equal(va, vb), f.name
        assert np.array_equal(a.h0.deaths, b.h0.deaths) and a.h0.n_essential == b.h0.n_essential
