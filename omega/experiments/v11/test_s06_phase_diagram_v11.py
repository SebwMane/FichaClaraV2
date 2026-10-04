"""Paso 6 / Experimento 3 (Ω-1.1): diagrama de fases reducido con taxonomia y certificado por punto (DESIGN §4).

COMPUERTA: en modo full, `require_prerequisites(out_root, "s06_phase_diagram")` exige O-00..O-04 en modo full completos
y se evalua ANTES de ejecutar nada. El modo smoke (N=24, 2x2 puntos) no tiene compuerta.

Preregistrado (R8, riesgo de expectativa honesta; escrito antes de medir y sin ajustar umbrales, M§44):
  * la dinamica S0 con clip da F0 para α̂ <= 1.9 y F1 para α̂ >= 2.1, para todo γ̂ de la malla (en el 100% de las corridas);
  * ningun punto es Ω-CANDIDATE y ninguna corrida pasa todos los campos;
  * la varianza entre semillas ("seed_variance", no termica) del peso medio es 0 en todo punto, de modo que λ* no
    esta definido (susceptibilidad maxima 0) para ningun γ̂.
Esta corrida usa un solo tamano: `size_robust` no es evaluable (F6 por evidencia faltante); el tamano finito es O-06.
Los nulos se evaluan solo en puntos con alguna corrida no trivial (primario fuera de F0/F1/F10).
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.certificate import point_verdict
from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.config.seeds import seed_key
from omega.config.settings11 import NullModel, Omega11Config
from omega.contracts import RunEvidence
from omega.controls.battery import null_battery, null_seed_key
from omega.experiments.v11.gate import (
    STEP_ORDER,
    PassportWriter,
    PrerequisiteError,
    assessment_row,
    derive_cfg,
    require_prerequisites,
    runs_root,
    write_summary,
)
from omega.network.weights import upper_triangle
from omega.phases.scan import default_config, point_params, simulate
from omega.statistics.transition import transition_summary

NAME = "s06_phase_diagram"
EXPERIMENT_ID = 1106
ENTRYPOINT = "omega.experiments.v11.test_s06_phase_diagram_v11:run"
ALPHA_SPLIT = 1.95
DEFAULT_ALPHAS = (0.0, 0.5, 1.0, 1.5, 1.9, 2.1, 2.5, 3.0, 4.0)
DEFAULT_GAMMAS = (0.0, 0.3, 1.0, 3.0, 10.0)
TRIVIAL = {"Ω-F0", "Ω-F1", "Ω-F10"}


def _grid(cfg: Omega11Config, mode: Literal["smoke", "full"]) -> tuple[int, int, tuple[float, ...], tuple[float, ...]]:
    if mode == "smoke":
        return 24, 2, (1.5, 3.0), (0.0, 10.0)
    scan = cfg.base.scan
    alphas = scan.alpha_hat if scan is not None else DEFAULT_ALPHAS
    gammas = scan.gamma_hat if scan is not None else DEFAULT_GAMMAS
    return 200, 10, tuple(alphas), tuple(gammas)


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta s06 y devuelve `summary.json`. En full exige O-00..O-04 full completos (PrerequisiteError si no)."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    n, reps, alphas, gammas = _grid(cfg, mode)
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=reps)
    expected = {
        "s0_codes": "F0 para α̂ <= 1.9 y F1 para α̂ >= 2.1 en el 100% de las corridas, para todo γ̂",
        "n_candidate_points": 0,
        "n_runs_pass": 0,
        "seed_variance_chi_max": 0.0,
        "lambda_star": "indefinido (chi_max == 0)",
    }
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID, "description": "diagrama de fases reducido, taxonomia y certificado por punto",
        "n": n, "replicates": reps, "alpha_hat": alphas, "gamma_hat": gammas, "expected": expected,
        "gate": "O-00..O-04 en modo full" if mode == "full" else "sin compuerta (smoke)",
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    w_min = cfg.base.graph.w_min
    points: list[dict[str, Any]] = []
    mean_weights: dict[float, dict[float, list[float]]] = {g: {a: [] for a in alphas} for g in gammas}
    idx = 0
    for a in alphas:
        for g in gammas:
            p = point_params(a, g, n)
            runs: list[RunEvidence] = []
            rows: list[dict[str, Any]] = []
            states = []
            for rep in range(reps):
                key = seed_key(cfg.base.seeds, idx, rep)
                res = simulate(cfg.base, p, key)
                w = res.trajectory.w_final
                c_ev = evidence_cfg(cfg, w)  # auditoria B1
                ev = collect_run_evidence(w, res.trajectory.status, c_ev, evidence_rng(key))
                runs.append(ev)
                states.append((key, w))
                row = assessment_row(ev, c_ev)
                row["omega10_label"] = res.assessment.label.value
                if rep == 0:
                    row["passport"] = writer.save(
                        derive_cfg(cfg, functional=p), seed=key, w=w, termination=res.trajectory.status.value,
                        init_distribution="uniform/upper_mirror", results=row, w0=res.w0)
                rows.append(row)
                mean_weights[g][a].append(float(upper_triangle(w).mean()))
            assessments = [assess_run(r, cfg) for r in runs]  # umbrales identicos en c_ev y cfg
            nontrivial = any(x.primary is not None and x.primary.value not in TRIVIAL for x in assessments) or any(
                x.passes for x in assessments)
            nulls: dict[NullModel, list[RunEvidence]] = {}
            if nontrivial:
                for rep, (key, w) in enumerate(states[: min(3, reps)]):
                    for kind, wn in null_battery(w, w_min, cfg.null_models, key).items():
                        nulls.setdefault(kind, []).append(
                            collect_run_evidence(wn, runs[rep].status, evidence_cfg(cfg, wn), evidence_rng(null_seed_key(key, kind))))
            pv = point_verdict(runs, {n: runs}, nulls, cfg)
            claim = "Ω-F0" if a < ALPHA_SPLIT else "Ω-F1"
            code_counts: Counter[str] = Counter(c for r in rows for c in r["codes"])
            points.append({
                "alpha_hat": a, "gamma_hat": g,
                "run_code_fractions": {k: v / reps for k, v in sorted(code_counts.items())},
                "run_outcome_fractions": dict(pv.run_outcome_fractions),
                "point_codes": [c.value for c in pv.codes],
                "point_primary": None if pv.primary is None else pv.primary.value,
                "verdict": pv.verdict.value,
                "failed_fields": list(pv.certificate.failed_fields),
                "consensus_dimension": pv.certificate.consensus_dimension,
                "dimension_class": pv.certificate.dimension_class,
                "certificate_evidence": dict(pv.certificate.evidence),
                "nulls_evaluated": bool(nulls),
                "claim_code": claim, "claim_met": code_counts[claim] == reps,
                "n_runs_pass": sum(1 for x in assessments if x.passes),
                "rows": rows,
            })
            idx += 1

    n_edges = n * (n - 1) // 2
    transition: dict[str, Any] = {}
    for g in gammas:
        samples = {a: np.asarray(v, dtype=np.float64) for a, v in mean_weights[g].items()}
        if len(samples) >= 2 and reps >= 2:
            ts = transition_summary(samples, n_edges)
            chi_max = float(np.max(ts.susceptibility))
            transition[f"{g:g}"] = {
                "label": "seed_variance", "chi_max": chi_max, "lambda_star": float(ts.peak_lambda) if chi_max > 0.0 else None,
                "lambdas": ts.lambdas, "mean": ts.mean, "susceptibility": ts.susceptibility,
            }
    candidates = [pt for pt in points if pt["verdict"] == "Ω-CANDIDATE"]
    summary = {
        **header, "stage": "final", "points": points, "seed_variance_transition": transition,
        "n_candidate_points": len(candidates), "n_runs_pass": sum(pt["n_runs_pass"] for pt in points),
        "s0_claims_met": all(pt["claim_met"] for pt in points),
        "lambda_star_undefined_everywhere": all(t["lambda_star"] is None for t in transition.values()),
        "expectations_met": (
            all(pt["claim_met"] for pt in points) and not candidates
            and all(pt["n_runs_pass"] == 0 for pt in points)
            and all(t["chi_max"] == 0.0 for t in transition.values())
        ),
        "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    assert len(data["points"]) == 4 and data["n_candidate_points"] == 0
    assert all(pt["verdict"] == "NOT_CANDIDATE" for pt in data["points"])
    assert set(data["seed_variance_transition"]) == {"0", "10"}
    assert len(data["passports"]) == 4
    assert (tmp_path / NAME / "summary.json").is_file()


def test_full_mode_requires_prerequisites(tmp_path: Path) -> None:
    """El modo full no ejecuta nada sin O-00..O-04 full (la compuerta se evalua primero)."""
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    with pytest.raises(PrerequisiteError):
        run(cfg, tmp_path, mode="full")
    assert not (tmp_path / NAME).exists()


def test_gate_semantics(tmp_path: Path) -> None:
    assert STEP_ORDER.index("s06_phase_diagram") == 5
    for step in STEP_ORDER[:5]:
        require_prerequisites(tmp_path, "o00_validation")  # el primero no exige nada
        write_summary(tmp_path, step, {"complete": True}, mode="smoke")
    with pytest.raises(PrerequisiteError):  # smoke no cuenta como full
        require_prerequisites(tmp_path, NAME)
    for step in STEP_ORDER[:5]:
        write_summary(tmp_path, step, {"complete": True}, mode="full")
    require_prerequisites(tmp_path, NAME)
    write_summary(tmp_path, "o03_nulls", {"complete": False}, mode="full")
    with pytest.raises(PrerequisiteError):
        require_prerequisites(tmp_path, NAME)
    with pytest.raises(ValueError):
        require_prerequisites(tmp_path, "paso_inexistente")


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=10))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita antes de ejecutar el paso 6
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
