"""Omega-8 (Ω-1.1, paso 9): topologia en estados no binarios y controles (DESIGN §1.2, §4).

Por familia (catalogo `state_families` de O-07): curvas de Betti(θ) (β0, β1^(4), densidad de β1, fraccion gigante), barras de
H0 (nacimientos y muertes de la filtracion), β1^(4) en W>w_min, complejo de cliques con χ (Euler por Betti y por conteos) y
`topology_stable`. Solo medida: la decision es del certificado.

Preregistrado (sin ajustar umbrales, M§44):
  * identidad de Euler: para todo complejo de cliques con status 'ok' y SIN simplices de dimension max_dim+1 (counts[-1] == 0, es
    decir, complejo de dimension <= 2 completamente enumerado), chi(Betti) == chi(conteos) (identidad matematica). Con 3-simplices
    el calculo trunca los Betti en dimension 2 y los conteos en 3, de modo que la igualdad no aplica (smoke lo mostro con ER/RGG;
    precondicion matematica anadida al preregistro tras el primer smoke, sin tocar ningun umbral);
  * arbol 3-ario: beta1^(4) == 0 (aciclico);
  * full (N~800): `topology_stable` True en anillo, T², T³ y RGG3 k12 (ciclos cortos de densidad <= 0.03) y False en ER k12
    (beta1^(4)/E ~ 0.3 >> 0.03: el cierre de 4-ciclos no rellena el espacio de ciclos aleatorio).
En smoke (N=40) las dos ultimas se informan pero no se afirman.
O-08 es el paso 9 del STEP_ORDER: en full exige O-00..O-05 y O-06 full completos.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.config.settings11 import Omega11Config
from omega.contracts import RunEvidence
from omega.experiments.v11.gate import (
    PassportWriter,
    assessment_row,
    check_expectation,
    derive_cfg,
    require_prerequisites,
    runs_root,
    write_summary,
)
from omega.experiments.v11.test_o07_dimensions import state_families
from omega.phases.finite_size import evidence_cfg
from omega.phases.scan import default_config
from omega.types import RunStatus

NAME = "o08_topology"
EXPERIMENT_ID = 1108
ENTRYPOINT = "omega.experiments.v11.test_o08_topology:run"
STABLE_TRUE = ("ring", "torus2d", "torus3d", "rgg3_k12")


def topology_row(ev: RunEvidence) -> dict[str, Any]:
    t = ev.topo
    c = t.curves
    deaths = t.h0.deaths
    clq = t.clique_at_w_min
    return {
        "betti_curves": {"thetas": c.thetas, "beta0": c.beta0, "beta1_short": c.beta1_short, "b1_density": c.b1_density,
                         "giant_fraction": c.giant_fraction, "status": list(c.status)},
        "h0": {"n_bars": int(t.h0.births.shape[0]), "n_essential": t.h0.n_essential,
               "death_min": float(deaths.min()) if deaths.size else None, "death_mean": float(deaths.mean()) if deaths.size else None,
               "death_max": float(deaths.max()) if deaths.size else None, "births": t.h0.births, "deaths": t.h0.deaths},
        "short_cycles_w_min": {"b0": t.at_w_min.b0, "b1": t.at_w_min.b1, "n_edges": t.at_w_min.n_edges, "n_faces": t.at_w_min.n_faces,
                               "b1_density": t.at_w_min.b1_density, "status": t.at_w_min.status},
        "clique_complex": {"betti": list(clq.betti), "counts": list(clq.counts), "euler_betti": clq.euler_betti,
                           "euler_counts": clq.euler_counts, "status": clq.status, "field": clq.field},
        "stable_flags": list(t.stable_flags), "topology_stable": bool(t.stable),
    }


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-08 y devuelve `summary.json`. En full exige los pasos 1..8 del STEP_ORDER (paso 9)."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    n = 40 if mode == "smoke" else 800
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=1)
    expected: dict[str, Any] = {
        "euler_identity_all_ok_complexes": True, "tree_b1_short": 0, "stable_true": list(STABLE_TRUE), "stable_false": ["er_k12"],
        "asserted": {"euler_identity": True, "tree_b1": True, "stable": mode == "full"},
    }
    header: dict[str, Any] = {"experiment": EXPERIMENT_ID, "description": "curvas de Betti, H0, beta1^(4), cliques y estabilidad topologica",
                              "n": n, "expected": expected}
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    for fam in state_families(cfg, n, mode):
        w, c_case, key = fam["w"], fam["cfg"], fam["key"]
        c_ev = evidence_cfg(c_case, w)
        st = RunStatus(fam["termination"]) if fam["termination"] in {s.value for s in RunStatus} else RunStatus.CONVERGED
        ev = collect_run_evidence(w, st, c_ev, evidence_rng(key))
        row = assessment_row(ev, c_ev)
        results[fam["name"]] = {"kind": fam["kind"], "n": int(w.shape[0]), **topology_row(ev), "assessment": row,
                                "passport": writer.save(c_case, seed=key, w=w, termination=fam["termination"],
                                                        init_distribution=fam["init"], results=row, null_model=fam["null_model"])}
    met: dict[str, bool] = {}
    ok_cliques = [r for r in results.values() if r["clique_complex"]["status"] == "ok" and r["clique_complex"]["counts"][-1] == 0]
    met["euler_identity"] = all(r["clique_complex"]["euler_betti"] == r["clique_complex"]["euler_counts"] for r in ok_cliques) and bool(ok_cliques)
    tree = next(k for k in results if k.startswith("tree3"))
    met["tree_b1_short"] = check_expectation("eq", results[tree]["short_cycles_w_min"]["b1"], 0)
    stable_met = all(check_expectation("eq", results[k]["topology_stable"], True) for k in STABLE_TRUE) and check_expectation(
        "eq", results["er_k12"]["topology_stable"], False)
    met["stable_flags_controls"] = stable_met
    asserted = ["euler_identity", "tree_b1_short"] + (["stable_flags_controls"] if mode == "full" else [])
    summary = {**header, "stage": "final", "results": results, "expectation_results": met,
               "expectations_met": all(met[k] for k in asserted), "passports": writer.labels, "complete": True}
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(40, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert {"ring", "torus2d", "torus3d", "rgg3_k12", "er_k12", "omega_u01_static"} <= set(res)
    assert data["expectation_results"]["euler_identity"] is True and data["expectation_results"]["tree_b1_short"] is True
    assert data["expectations_met"] is True
    r = res["ring"]
    assert r["short_cycles_w_min"]["b1"] == 1 and r["h0"]["n_essential"] == 1 and len(r["betti_curves"]["thetas"]) == 11
    assert len(data["passports"]) == len(res)
    assert np.isfinite(res["torus2d"]["short_cycles_w_min"]["b1_density"])


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=1))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita: paso 9 del STEP_ORDER
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
