"""Omega-10 (Ω-1.1, paso 9): curvatura de Ollivier (media y cola) y perfil QRC, y control S² frente a T² (DESIGN §1.6, §4).

Medidas por familia (catalogo `state_families` de O-07, a las que se anade la evidencia con QRC): kappa medio, desviacion, SE,
fraccion de cola kappa < -0.5, aristas muestreadas y `curvature_ok`; el perfil QRC (d̄(S_δ(x),S_δ(y))/δ, δ=1,2,3) es SOLO
INFORME y no entra en el certificado.

Preregistrado (sin ajustar umbrales, M§44):
  * controles exactos (§1.6): anillo, T² y T³ hipercubicos con lados >= 6 tienen kappa = 0 en cada arista (|kappa medio| <= 1e-9,
    cola 0). Se afirma tambien en smoke;
  * control S² frente a T² repetido (RGG2 binario, k=10, mismas reps): <kappa>(S²) > <kappa>(T²) (curvatura positiva frente a plana;
    la senal es ~1% de la escala y exige N grande, R1). Se afirma solo en full (N=800, 5 reps); en smoke (N=100, 2 reps y 50 aristas)
    se informa.
O-10 es el paso 10 del STEP_ORDER (indice 9): en full exige O-00..O-08 full completos.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.config.seeds import make_rng, seed_key
from omega.config.settings11 import Omega11Config
from omega.contracts import RunEvidence
from omega.controls.random_geometric import rgg_sphere, rgg_torus
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
from omega.experiments.v11.test_o07_dimensions import state_families
from omega.phases.finite_size import evidence_cfg, size_config
from omega.phases.scan import default_config
from omega.types import RunStatus

NAME = "o10_curvature"
EXPERIMENT_ID = 1110
ENTRYPOINT = "omega.experiments.v11.test_o10_curvature:run"
EXACT_TOL = 1e-9


def curvature_row(ev: RunEvidence) -> dict[str, Any]:
    c = ev.curvature
    row: dict[str, Any] = {"kappa_mean": c.mean, "kappa_std": c.std, "kappa_se": c.se, "tail_fraction": c.tail_fraction,
                           "n_edges": c.n_edges, "sampled": c.sampled, "curvature_ok": c.ok, "qrc": None}
    if ev.qrc is not None:
        q = ev.qrc
        row["qrc"] = {"deltas": q.deltas, "ratio_mean": q.ratio_mean, "ratio_se": q.ratio_se, "n_pairs": q.n_pairs,
                      "note": "solo informe; no entra en el certificado"}
    return row


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-10 y devuelve `summary.json`. En full exige los pasos 1..9 del STEP_ORDER."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    n, n_sph, reps, max_edges = (50, 100, 2, 50) if mode == "smoke" else (800, 800, 5, cfg.curvature.max_edges)
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=max(reps, 1))
    cfg = dataclasses.replace(cfg, curvature=dataclasses.replace(cfg.curvature, max_edges=max_edges))
    exact_shapes: dict[str, tuple[int, ...]] = (
        {"ring": (24,), "torus2d": (6, 6), "torus3d": (6, 6, 6)} if mode == "smoke"
        else {"ring": (800,), "torus2d": (28, 28), "torus3d": (9, 9, 9)}
    )
    expected: dict[str, Any] = {
        "exact_kappa_zero": list(exact_shapes), "exact_tol": EXACT_TOL,
        "sphere_gt_torus": "<kappa>(S²) > <kappa>(T²)", "asserted_sphere_gt_torus": mode == "full",
        "qrc": "solo informe",
    }
    header: dict[str, Any] = {"experiment": EXPERIMENT_ID, "description": "curvatura de Ollivier, QRC (informe) y control S² frente a T²",
                              "n": n, "max_edges": max_edges, "expected": expected}
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode, cfg=cfg)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    for fam in state_families(cfg, n, mode):
        w, c_case, key = fam["w"], size_config(fam["cfg"], fam["w"].shape[0]), fam["key"]
        c_ev = evidence_cfg(dataclasses.replace(c_case, curvature=cfg.curvature), w)
        st = RunStatus(fam["termination"]) if fam["termination"] in {s.value for s in RunStatus} else RunStatus.CONVERGED
        ev = collect_run_evidence(w, st, c_ev, evidence_rng(key), with_qrc=True)
        row = assessment_row(ev, c_ev)
        results[fam["name"]] = {"kind": fam["kind"], "n": int(w.shape[0]), **curvature_row(ev), "assessment": row,
                                "passport": writer.save(c_case, seed=key, w=w, termination=fam["termination"],
                                                        init_distribution=fam["init"], results=row, null_model=fam["null_model"])}

    exact: dict[str, Any] = {}
    for i, (name, shape) in enumerate(exact_shapes.items()):
        w = periodic_lattice(shape)
        key = seed_key(cfg.base.seeds, 100 + i, 0)
        c_case = dataclasses.replace(size_config(cfg, int(w.shape[0])), curvature=cfg.curvature)
        ev = collect_run_evidence(w, RunStatus.CONVERGED, c_case, evidence_rng(key))
        exact[name] = {"shape": list(shape), **curvature_row(ev)}

    pairs: list[dict[str, Any]] = []
    for rep in range(reps):
        vals: dict[str, float] = {}
        for tag in ("sphere", "torus"):
            key = seed_key(cfg.base.seeds, 200 + (0 if tag == "sphere" else 1), rep)
            rng = make_rng(key)
            ws = rgg_sphere(n_sph, 2, 10.0, rng) if tag == "sphere" else rgg_torus(n_sph, 2, 10.0, rng)
            c_case = dataclasses.replace(size_config(cfg, n_sph), curvature=cfg.curvature)
            ev = collect_run_evidence(ws, RunStatus.CONVERGED, c_case, evidence_rng(key))
            vals[tag] = ev.curvature.mean
            if rep == 0:
                writer.save(c_case, seed=key, w=ws, termination="STATIC", init_distribution=f"control/rgg2_{tag}",
                            results=curvature_row(ev))
        pairs.append({"rep": rep, "kappa_sphere": vals["sphere"], "kappa_torus": vals["torus"], "diff": vals["sphere"] - vals["torus"]})
    diffs = np.asarray([p["diff"] for p in pairs], dtype=np.float64)
    sphere_test = {"mean_diff": float(diffs.mean()), "se": float(diffs.std(ddof=1) / np.sqrt(diffs.size)) if diffs.size > 1 else None,
                   "n": int(diffs.size), "pairs": pairs, "n_sph": n_sph}

    met = {f"{k}.kappa_zero": bool(check_expectation("abs_within", v["kappa_mean"], 0.0, EXACT_TOL)
                                   and check_expectation("eq", v["tail_fraction"], 0.0)) for k, v in exact.items()}
    met["sphere_gt_torus"] = check_expectation("gt", sphere_test["mean_diff"], 0.0)
    asserted = list(k for k in met if k != "sphere_gt_torus") + (["sphere_gt_torus"] if mode == "full" else [])
    summary = {**header, "stage": "final", "results": results, "exact_controls": exact, "sphere_vs_torus": sphere_test,
               "expectation_results": met, "expectations_met": all(met[k] for k in asserted),
               "passports": writer.labels, "complete": True}
    return write_summary(out_root, NAME, summary, mode=mode, cfg=cfg)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(50, master_entropy=20240901, replicates=2))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    ex = data["exact_controls"]
    assert set(ex) == {"ring", "torus2d", "torus3d"}
    for k, v in ex.items():
        assert abs(v["kappa_mean"]) <= EXACT_TOL and v["tail_fraction"] == 0.0, k
    assert data["expectation_results"]["ring.kappa_zero"] is True and data["expectations_met"] is True
    assert data["sphere_vs_torus"]["n"] == 2
    assert all(r["n_edges"] <= 50 for r in data["results"].values())
    assert data["results"]["ring"]["qrc"] is not None and len(data["passports"]) >= len(data["results"])


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=5))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita: paso 10 del STEP_ORDER
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
