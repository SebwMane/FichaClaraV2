"""Omega-4 (Ω-1.1): ablacion del funcional y rama Omega-B de densidad fija (DESIGN §1.9, §4, riesgo R8).

Preregistrado:
  * NO_TRIANGLES (T=0, sin termino de triangulos) da F0 en el 100% de las corridas (analitico: el gradiente es >= 0);
  * TRIANGLES_ONLY da F1 en el 100% (el gradiente es <= 0 y el flujo llega a J-I);
  * Omega-B, α̂ = 0.5·α̂_c(γ̂, ρ): estado uniforme (F1) en el 100%;
  * Omega-B, γ̂=0 y α̂ >= 1.5·α̂_c: clique mas aislados, es decir F2 en el 100%.
  El resto de celdas (FULL, NO_DENSITY, NO_DEGREE, Omega-B con γ̂>0 sobre el umbral) es solo informe.
"El codigo X" significa X en el conjunto completo de codigos de la corrida.
Full: N=100 (no se ejecuta N=200 por coste O(N^3) por paso), 10 semillas. Smoke: N=24.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Literal

import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.config.convert import reduced_to_raw
from omega.config.seeds import make_rng, seed_key
from omega.config.settings11 import Engine, FixedDensityConfig, Omega11Config
from omega.dynamics.fixed_density import evolve_fixed_density, uniform_state_threshold
from omega.experiments.v11.gate import PassportWriter, assessment_row, derive_cfg, require_prerequisites, runs_root, write_summary
from omega.network.initialization import random_uniform_weights
from omega.phases.ablation import AblationSpec, ablation_set, evolve_ablated
from omega.phases.finite_size import evidence_cfg
from omega.phases.scan import default_config

NAME = "o04_ablation"
EXPERIMENT_ID = 1104
ENTRYPOINT = "omega.experiments.v11.test_o04_ablation:run"

ABLATION_CLAIM = {AblationSpec.NO_TRIANGLES: "Ω-F0", AblationSpec.TRIANGLES_ONLY: "Ω-F1"}
EXTRA_SPECS = (AblationSpec.FULL, AblationSpec.NO_TRIANGLES)  # Enmienda A-2: tambien con el N mayor (auditoria, desviacion 7)


def _grid(mode: Literal["smoke", "full"]) -> dict[str, Any]:
    if mode == "smoke":
        return {"n": 24, "reps": 1, "alphas": (0.5, 4.0), "gammas": (0.0, 10.0),
                "rhos": (0.1,), "factors": (0.5, 3.0), "b_gammas": (0.0, 1.0), "n_extra": 28}
    return {"n": 100, "reps": 10, "alphas": (0.5, 1.5, 2.5, 4.0), "gammas": (0.0, 1.0, 10.0),
            "rhos": (0.05, 0.1, 0.2), "factors": (0.5, 1.5, 3.0), "b_gammas": (0.0, 1.0, 10.0), "n_extra": 200}


def _claim_b(factor: float, gamma: float) -> str | None:
    if factor < 1.0:
        return "Ω-F1"
    return "Ω-F2" if gamma == 0.0 else None


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-04 y devuelve `summary.json`. En full exige O-00..O-03 full completos."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    g = _grid(mode)
    n, reps = int(g["n"]), int(g["reps"])
    sizes = (n, int(g["n_extra"]))
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=reps)
    expected: list[dict[str, Any]] = []
    for spec, code in ABLATION_CLAIM.items():
        for n_i in sizes if spec in EXTRA_SPECS else (n,):
            expected.append({"block": "ablation", "spec": spec.value, "n": n_i, "code": code, "fraction": 1.0,
                             "alpha_hat": list(g["alphas"]), "gamma_hat": list(g["gammas"])})
    for rho in g["rhos"]:
        for gm in g["b_gammas"]:
            for f in g["factors"]:
                c = _claim_b(f, gm)
                if c:
                    expected.append({"block": "omega_b", "rho": rho, "gamma_hat": gm, "factor": f, "code": c, "fraction": 1.0})
    header: dict[str, Any] = {"experiment": EXPERIMENT_ID, "description": "ablacion y Omega-B (densidad fija)",
                              "n": n, "n_extra": sizes[1], "extra_specs": [x.value for x in EXTRA_SPECS],
                              "replicates": reps, "expected": expected}
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    idx = 0
    cells: list[dict[str, Any]] = []
    for a in g["alphas"]:
        for gm in g["gammas"]:
            for n_i in sizes:
                c_n = derive_cfg(cfg, n=n_i)
                full_p = reduced_to_raw(a, gm, n_i)
                for ap in ablation_set(a, gm, n_i):
                    if n_i != n and ap.spec not in EXTRA_SPECS:
                        continue
                    rows: list[dict[str, Any]] = []
                    for rep in range(reps):
                        key = seed_key(cfg.base.seeds, idx, rep)
                        w0 = random_uniform_weights(n_i, make_rng(key))
                        traj = evolve_ablated(w0, ap, cfg.base.dynamics)
                        c_ev = evidence_cfg(c_n, traj.w_final)  # auditoria B1
                        ev = collect_run_evidence(traj.w_final, traj.status, c_ev, evidence_rng(key))
                        row = assessment_row(ev, c_ev)
                        row["steps"] = traj.steps
                        row["passport"] = writer.save(
                            derive_cfg(c_n, functional=full_p), seed=key, w=traj.w_final, termination=traj.status.value,
                            init_distribution=f"uniform/upper_mirror;ablation={ap.spec.value}",
                            results={**row, "ablation": ap.spec.value}, w0=w0)  # auditoria B8/B15
                        rows.append(row)
                    idx += 1
                    counts: Counter[str] = Counter(c for r in rows for c in r["codes"])
                    claim = ABLATION_CLAIM.get(ap.spec)
                    cells.append({
                        "block": "ablation", "spec": ap.spec.value, "n": n_i, "alpha_hat": a, "gamma_hat": gm,
                        "code_fractions": {k: v / reps for k, v in sorted(counts.items())},
                        "primary_fractions": {k: v / reps for k, v in sorted(Counter(str(r["primary"]) for r in rows).items())},
                        "claim_code": claim, "claim_met": None if claim is None else counts[claim] == reps,
                        "n_passes": sum(1 for r in rows if r["passes"]), "rows": rows,
                    })

    for rho in g["rhos"]:
        for gm in g["b_gammas"]:
            alpha_c = uniform_state_threshold(gm, rho, n)
            for f in g["factors"]:
                a = f * alpha_c
                p = reduced_to_raw(a, gm, n)
                fd = FixedDensityConfig(rho=rho)
                c_case = derive_cfg(cfg, functional=p, engine=Engine.FIXED_DENSITY, fixed_density=fd)
                rows = []
                for rep in range(reps):
                    key = seed_key(cfg.base.seeds, idx, rep)
                    w0 = random_uniform_weights(n, make_rng(key))
                    traj = evolve_fixed_density(w0, p, cfg.base.dynamics, fd)
                    c_ev = evidence_cfg(c_case, traj.w_final)  # auditoria B1
                    ev = collect_run_evidence(traj.w_final, traj.status, c_ev, evidence_rng(key))
                    row = assessment_row(ev, c_ev)
                    row["steps"] = traj.steps
                    row["passport"] = writer.save(c_case, seed=key, w=traj.w_final, termination=traj.status.value,
                                                  init_distribution="uniform/upper_mirror;fixed_density", results=row, w0=w0)
                    rows.append(row)
                idx += 1
                counts = Counter(c for r in rows for c in r["codes"])
                claim_b = _claim_b(f, gm)
                cells.append({
                    "block": "omega_b", "rho": rho, "gamma_hat": gm, "factor": f, "alpha_hat": a, "alpha_c": alpha_c,
                    "code_fractions": {k: v / reps for k, v in sorted(counts.items())},
                    "primary_fractions": {k: v / reps for k, v in sorted(Counter(str(r["primary"]) for r in rows).items())},
                    "claim_code": claim_b, "claim_met": None if claim_b is None else counts[claim_b] == reps,
                    "n_passes": sum(1 for r in rows if r["passes"]), "rows": rows,
                })
    judged = [c["claim_met"] for c in cells if c["claim_met"] is not None]
    summary = {
        **header, "stage": "final", "cells": cells,
        "n_claims": len(judged), "n_claims_met": sum(1 for j in judged if j), "expectations_met": all(judged),
        "n_runs_pass": sum(c["n_passes"] for c in cells), "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    cells = data["cells"]
    assert sum(1 for c in cells if c["block"] == "ablation") == (5 + 2) * 2 * 2  # N base + FULL/NO_TRIANGLES en n_extra
    assert {c["n"] for c in cells if c["block"] == "ablation" and c["spec"] == "full"} == {24, 28}
    assert {c["n"] for c in cells if c["block"] == "ablation" and c["spec"] == "no_density"} == {24}
    assert sum(1 for c in cells if c["block"] == "omega_b") == 1 * 2 * 2
    # T=0 da F0 siempre (analitico) y TRIANGLES_ONLY da F1.
    for c in cells:
        if c["block"] == "ablation" and c["spec"] == "no_triangles":
            assert c["claim_met"] is True, c
        if c["block"] == "ablation" and c["spec"] == "triangles_only":
            assert c["claim_met"] is True, c
    assert len(data["passports"]) == len(cells) * data["replicates"]  # un pasaporte por corrida (B8)
    ab = [json.loads(p.read_text(encoding="utf-8")) for p in (tmp_path / NAME / "passports").glob("OMEGA-EXP-*.json")]
    assert len(ab) == len(data["passports"])
    assert {p["results"].get("ablation") for p in ab if "ablation" in p["results"]} == {x.value for x in AblationSpec}  # B15


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(100, master_entropy=20240901, replicates=10))
    data = json.loads(run(cfg, runs_root(), mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
