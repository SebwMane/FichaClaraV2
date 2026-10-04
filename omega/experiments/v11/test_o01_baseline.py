"""Omega-1 (Ω-1.1): linea base S0 de Omega-1.0 reevaluada con la taxonomia v1.1 (DESIGN §4, riesgo R8).

Preregistrado: con N=200, α̂ ∈ {0.5, 1, 1.5, 1.9, 2.1, 3}, γ̂ ∈ {0, 1, 10} y 10 semillas, α̂ <= 1.9 da F0 en el 100% de las
corridas y α̂ >= 2.1 da F1 en el 100% (registro negativo: S0 no produce regiones geometricas). "Da F0/F1" significa que el
codigo pertenece al conjunto completo de codigos de la corrida; el primario se informa aparte (F10 lo precede si no converge).
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Literal

import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.config.seeds import seed_key
from omega.config.settings11 import Omega11Config
from omega.experiments.v11.gate import PassportWriter, assessment_row, derive_cfg, require_prerequisites, runs_root, write_summary
from omega.phases.finite_size import evidence_cfg
from omega.phases.scan import default_config, point_params, simulate

NAME = "o01_baseline"
EXPERIMENT_ID = 1101
ENTRYPOINT = "omega.experiments.v11.test_o01_baseline:run"
ALPHA_SPLIT = 1.95  # α̂ <= 1.9 (F0) frente a α̂ >= 2.1 (F1)


def _grid(mode: Literal["smoke", "full"]) -> tuple[int, int, tuple[float, ...], tuple[float, ...]]:
    if mode == "smoke":
        return 32, 2, (1.5, 3.0), (0.0, 10.0)
    return 200, 10, (0.5, 1.0, 1.5, 1.9, 2.1, 3.0), (0.0, 1.0, 10.0)


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-01 y devuelve `summary.json`. En full exige O-00 full completo."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    n, reps, alphas, gammas = _grid(mode)
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=reps)
    expected = [
        {"id": f"a{a:g}_g{g:g}", "alpha_hat": a, "gamma_hat": g, "claim": "F0 en 100%" if a < ALPHA_SPLIT else "F1 en 100%",
         "code": "Ω-F0" if a < ALPHA_SPLIT else "Ω-F1", "fraction": 1.0}
        for a in alphas for g in gammas
    ]
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID,
        "description": "S0 de Omega-1.0 con taxonomia v1.1 (registro negativo, riesgo R8)",
        "n": n, "replicates": reps, "alpha_hat": alphas, "gamma_hat": gammas, "expected": expected,
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode, cfg=cfg)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    points: list[dict[str, Any]] = []
    idx = 0
    for a in alphas:
        for g in gammas:
            p = point_params(a, g, n)
            rows: list[dict[str, Any]] = []
            for rep in range(reps):
                key = seed_key(cfg.base.seeds, idx, rep)
                res = simulate(cfg.base, p, key)
                w = res.trajectory.w_final
                c_ev = evidence_cfg(cfg, w)  # auditoria B1
                ev = collect_run_evidence(w, res.trajectory.status, c_ev, evidence_rng(key))
                row = assessment_row(ev, c_ev)
                row["omega10_label"] = res.assessment.label.value
                row["steps"] = res.trajectory.steps
                rows.append(row)
                row["passport"] = writer.save(
                    derive_cfg(cfg, functional=p), seed=key, w=w, termination=res.trajectory.status.value,
                    init_distribution="uniform/upper_mirror", results=row, w0=res.w0)
            code_counts: Counter[str] = Counter(c for r in rows for c in r["codes"])
            claim = "Ω-F0" if a < ALPHA_SPLIT else "Ω-F1"
            points.append({
                "id": f"a{a:g}_g{g:g}", "alpha_hat": a, "gamma_hat": g,
                "code_fractions": {k: v / reps for k, v in sorted(code_counts.items())},
                "primary_fractions": {k: v / reps for k, v in sorted(Counter(str(r["primary"]) for r in rows).items())},
                "omega10_labels": dict(Counter(r["omega10_label"] for r in rows)),
                "claim_code": claim, "claim_met": code_counts[claim] == reps,
                "n_passes": sum(1 for r in rows if r["passes"]), "rows": rows,
            })
            idx += 1
    summary = {
        **header, "stage": "final", "points": points,
        "expectations_met": all(p["claim_met"] for p in points),
        "n_points_claim_met": sum(1 for p in points if p["claim_met"]),
        "n_runs_pass": sum(p["n_passes"] for p in points),
        "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode, cfg=cfg)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(32, master_entropy=20240901, replicates=2))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    assert len(data["points"]) == 4 and len(data["expected"]) == 4
    assert all(len(p["rows"]) == 2 for p in data["points"])
    assert data["n_runs_pass"] == 0  # S0 no produce candidatos (R8)
    assert len(data["passports"]) == 4 * 2  # un pasaporte por corrida dinamica (auditoria B8)
    assert len(list((tmp_path / NAME / "passports").glob("OMEGA-EXP-*.json"))) == 4 * 2


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=10))
    data = json.loads(run(cfg, runs_root(), mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True and len(data["points"]) == 18
