"""Omega-7 (Ω-1.1): comparacion de estimadores de dimension con curvas por escala (DESIGN §1.1, §1.3, §4).

Para cada familia (estados Ω y controles) se miden D_vol, D_ball, D_hop, D_log, D_inv, D_res (suite de distancias), D_s y D_Weyl
con sus curvas por escala (escalas, perfil y pendientes locales) y el acuerdo por pares |ΔD| <= `certificate.class_tol` (0.35).
O-07 no tiene orden propio en el STEP_ORDER (§4): su modo full exige lo mismo que O-08.

Este modulo tambien define `state_families`, el catalogo de estados que comparten O-07, O-08, O-09 y O-10:
  * Ω: S0 (α̂=1.5 y 3, γ̂=1; N_Ω = min(N, 200)), U(0,1) estatico y un estado Langevin (Θ̂=0.3, α̂=1.5, γ̂=1);
  * controles: anillo, T², T³, RGG2 k10, RGG3 k12, Watts-Strogatz p=0.05, ER k12 y arbol 3-ario.

Preregistrado (full, N~800; valores medidos en DESIGN §1.3-§1.4, sin ajustar umbrales, M§44):
  * clase de dimension del consenso: anillo 1, T² 2, T³ 3, RGG3 k12 3 y los tres estimadores D_vol, D_s, D_Weyl concuerdan por pares
    (<= 0.35) en T² y T³ y en RGG3;
  * ER k12 contiene el codigo Ω-F4 (dimension no definida/artefacto);
  * estados Ω (S0, U(0,1), Langevin): solo informe.
En smoke (N=40) las expectativas se evaluan y se informan pero no se afirman (la ventana de escala exige N grande, R1).
"""

from __future__ import annotations

import dataclasses
import itertools
import json
import math
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings11 import DistanceMode, LangevinConfig, NullModel, Omega11Config
from omega.contracts import FailureCode, RunEvidence
from omega.controls.random_geometric import balanced_tree, rgg_torus
from omega.controls.small_world import watts_strogatz
from omega.controls.erdos_renyi import erdos_renyi_gnp
from omega.experiments.reference_graphs import periodic_lattice
from omega.experiments.v11.gate import (
    PassportWriter,
    assessment_row,
    check_expectation,
    derive_cfg,
    require_prerequisites,
    runs_root,
    write_summary,
)
from omega.network.initialization import random_uniform_weights
from omega.phases.finite_size import evidence_cfg, size_config
from omega.phases.scan import default_config, point_params, simulate
from omega.statistics.langevin import run_langevin
from omega.types import DimensionEstimate, FloatArray, RunStatus

__all__ = ["state_families", "estimator_table", "NAME", "EXPERIMENT_ID"]

NAME = "o07_dimensions"
EXPERIMENT_ID = 1107
ENTRYPOINT = "omega.experiments.v11.test_o07_dimensions:run"
PREREQ_STEP = "o08_topology"  # O-07 y O-09 no tienen orden propio: requieren lo mismo que O-08


def state_families(cfg: Omega11Config, n: int, mode: Literal["smoke", "full"]) -> list[dict[str, Any]]:
    """Catalogo de estados (Ω y controles). Cada elemento: name, kind ('omega'|'control'), w, cfg (del caso), key,
    termination, null_model, init. `cfg.base.seeds.replicates` debe ser >= 1; la clave de la familia i es seed_key(.., i, 0)."""
    n_om = min(n, 200) if mode == "full" else min(n, 24)
    fams: list[dict[str, Any]] = []

    def add(name: str, kind: str, w: FloatArray, key: SeedKey, termination: str, init: str,
            null: NullModel | None = None, c: Omega11Config | None = None) -> None:
        base = c if c is not None else cfg
        fams.append({"name": name, "kind": kind, "w": w, "cfg": size_config(base, int(w.shape[0])), "key": key,
                     "termination": termination, "null_model": null, "init": init})

    i = 0

    def nk() -> SeedKey:
        nonlocal i
        k = seed_key(cfg.base.seeds, i, 0)
        i += 1
        return k

    for a in (1.5, 3.0):
        key = nk()
        p = point_params(a, 1.0, n_om)
        c0 = derive_cfg(size_config(cfg, n_om), functional=p)
        res = simulate(c0.base, p, key)
        add(f"omega_s0_a{a:g}_g1", "omega", res.trajectory.w_final, key, res.trajectory.status.value, "uniform/upper_mirror", c=c0)
    key = nk()
    add("omega_u01_static", "omega", random_uniform_weights(n_om, make_rng(key)), key, "STATIC", "uniform/upper_mirror")
    key = nk()
    p = point_params(1.5, 1.0, n_om)
    lc = LangevinConfig(theta_hat=0.3, n_steps=100 if mode == "smoke" else 2000)
    w0 = random_uniform_weights(n_om, make_rng(key))
    chain = run_langevin(w0, p, lc, cfg.base.graph.w_min, make_rng(SeedKey(key.entropy, key.spawn_key + (3000,))))
    add("omega_langevin_th0.3_a1.5_g1", "omega", chain.w_final, key, "CHAIN_COMPLETE", "uniform/upper_mirror;engine=langevin",
        c=derive_cfg(size_config(cfg, n_om), functional=p))

    s2 = max(4, round(math.sqrt(n)))
    s3 = max(4, round(n ** (1.0 / 3.0)))
    h = 3 if n < 200 else 5
    for name, w, null in (
        ("ring", periodic_lattice((n,)), None),
        ("torus2d", periodic_lattice((s2, s2)), None),
        ("torus3d", periodic_lattice((s3, s3, s3)), None),
    ):
        add(name, "control", w, nk(), "STATIC", f"control/{name}", null)
    key = nk()
    add("rgg2_k10", "control", rgg_torus(n, 2, 10.0, make_rng(key)), key, "STATIC", "control/rgg2_k10", NullModel.RANDOM_GEOMETRIC)
    key = nk()
    add("rgg3_k12", "control", rgg_torus(n, 3, 12.0, make_rng(key)), key, "STATIC", "control/rgg3_k12", NullModel.RANDOM_GEOMETRIC)
    key = nk()
    add("ws_p0.05", "control", watts_strogatz(n, 12, 0.05, make_rng(key)), key, "STATIC", "control/ws_p0.05", NullModel.SMALL_WORLD)
    key = nk()
    add("er_k12", "control", erdos_renyi_gnp(n, 12.0 / (n - 1), make_rng(key)), key, "STATIC", "control/er_k12", NullModel.ERDOS_RENYI)
    add(f"tree3_h{h}", "control", balanced_tree(3, h), nk(), "STATIC", f"control/tree3_h{h}")
    return fams


def _status(s: str) -> RunStatus:
    return RunStatus(s) if s in {x.value for x in RunStatus} else RunStatus.CONVERGED


def _est(e: DimensionEstimate) -> dict[str, Any]:
    return {"value": e.value, "stderr": e.stderr, "status": e.status, "window": e.window, "plateau": e.plateau, "method": e.method,
            "curve": {"scales": e.scales, "profile": e.profile, "local_slopes": e.local_slopes}}


def estimator_table(ev: RunEvidence) -> dict[str, DimensionEstimate]:
    """Estimadores de dimension de una corrida, por nombre."""
    geo, suite = ev.geometry, ev.suite
    table: dict[str, DimensionEstimate] = {"d_vol": geo.d_eff, "d_ball": geo.d_eff_ball, "d_hop": geo.d_eff_hops,
                                           "d_s": geo.d_s, "d_weyl": ev.d_weyl}
    names = {DistanceMode.WEIGHTED_LOG: "d_log", DistanceMode.WEIGHTED_INVERSE: "d_inv", DistanceMode.RESISTANCE: "d_res"}
    for mode, nm in names.items():
        if mode in suite.estimates:
            table[nm] = suite.estimates[mode]
    return table


def measure(ev: RunEvidence, cfg: Omega11Config) -> dict[str, Any]:
    table = estimator_table(ev)
    tol = cfg.certificate.class_tol
    pairs: dict[str, bool | None] = {}
    for (a, ea), (b, eb) in itertools.combinations(table.items(), 2):
        ok = math.isfinite(ea.value) and math.isfinite(eb.value)
        pairs[f"{a}|{b}"] = bool(abs(ea.value - eb.value) <= tol) if ok else None
    core = ("d_vol", "d_s", "d_weyl")
    core_ok = [pairs.get(f"{a}|{b}") for a, b in itertools.combinations(core, 2)]
    return {
        "estimators": {k: _est(v) for k, v in table.items()},
        "pairwise_agree": pairs, "core_pairwise_agree": all(x is True for x in core_ok),
        "consensus_dimension": ev.consensus_dimension, "dimension_class": ev.dimension_class,
        "metric_spread": ev.suite.metric_spread, "resistance_exponent": ev.suite.resistance_exponent,
    }


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-07 y devuelve `summary.json`. En full exige lo mismo que O-08 (sin orden propio)."""
    if mode == "full":
        require_prerequisites(out_root, PREREQ_STEP)
    n = 40 if mode == "smoke" else 800
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=1)
    expected: dict[str, Any] = {
        "class": {"ring": 1, "torus2d": 2, "torus3d": 3, "rgg3_k12": 3},
        "core_pairwise_agree": ["torus2d", "torus3d", "rgg3_k12"],
        "er_k12_contains": FailureCode.F4.value,
        "omega_states": "solo informe",
        "asserted": mode == "full",
    }
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID, "description": "estimadores de dimension y acuerdo por pares (curvas por escala)",
        "n": n, "expected": expected, "prerequisites": PREREQ_STEP if mode == "full" else "sin compuerta (smoke)",
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode, cfg=cfg)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    for fam in state_families(cfg, n, mode):
        w, c_case, key = fam["w"], fam["cfg"], fam["key"]
        c_ev = evidence_cfg(c_case, w)
        ev = collect_run_evidence(w, _status(fam["termination"]), c_ev, evidence_rng(key))
        row = assessment_row(ev, c_ev)
        results[fam["name"]] = {"kind": fam["kind"], "n": int(w.shape[0]), **measure(ev, c_ev), "assessment": row,
                                "passport": writer.save(c_case, seed=key, w=w, termination=fam["termination"],
                                                        init_distribution=fam["init"], results=row, null_model=fam["null_model"])}
    exp = expected
    met: dict[str, bool] = {}
    for name, cls in exp["class"].items():
        met[f"{name}.class"] = check_expectation("eq", results[name]["dimension_class"], cls)
    for name in exp["core_pairwise_agree"]:
        met[f"{name}.core_pairwise_agree"] = check_expectation("eq", results[name]["core_pairwise_agree"], True)
    er = next(k for k in results if k.startswith("er_k12"))
    met["er_k12.F4"] = check_expectation("contains", results[er]["assessment"]["codes"], FailureCode.F4.value)
    summary = {**header, "stage": "final", "results": results, "expectation_results": met,
               "expectations_met": all(met.values()), "passports": writer.labels, "complete": True}
    return write_summary(out_root, NAME, summary, mode=mode, cfg=cfg)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(40, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert {"ring", "torus2d", "torus3d", "rgg3_k12", "er_k12", "omega_u01_static", "omega_s0_a1.5_g1"} <= set(res)
    for r in res.values():
        assert {"d_vol", "d_ball", "d_hop", "d_s", "d_weyl"} <= set(r["estimators"])
        assert "curve" in r["estimators"]["d_vol"] and isinstance(r["core_pairwise_agree"], bool)
    assert len(data["passports"]) == len(res)
    assert isinstance(data["expectations_met"], bool)


def test_state_families_catalog() -> None:
    cfg = Omega11Config(base=default_config(40, master_entropy=20240901, replicates=1))
    fams = state_families(cfg, 40, "smoke")
    names = [f["name"] for f in fams]
    assert len(set(names)) == len(names) == 12
    assert {f["kind"] for f in fams} == {"omega", "control"}
    assert all(f["cfg"].base.init.n == f["w"].shape[0] for f in fams)
    again = state_families(cfg, 40, "smoke")
    assert all(np.array_equal(a["w"], b["w"]) for a, b in zip(fams, again, strict=True))  # determinista por semilla
    assert dataclasses.is_dataclass(fams[0]["cfg"])


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=1))
    out = runs_root()
    require_prerequisites(out, PREREQ_STEP)  # sin orden propio: lo mismo que O-08
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
