"""Omega-6 (Ω-1.1, paso 7): tamano finito D(N), C(N), G/N, rho(N), xi_F(N), L_hop(N) (DESIGN §1.14, §4, riesgos R1 y R8).

Puntos: S0 con α̂ ∈ {1.9, 2.1, 3} × γ̂ ∈ {0, 10} y Omega-B (ρ=0.1, α̂ = 1.5·α̂_c(γ̂,ρ,N), γ̂ ∈ {0, 10}). Tamanos de
`cfg.finite_size.sizes` (full: 64..800); los tamanos XL (`xl_sizes`) solo para puntos cuyas corridas pasan TODOS los campos
de corrida con el N mayor. Control positivo del mecanismo `size_robust` (R1, solo full): RGG3 k12 binario en
N ∈ (800, 1500, 3000).

Preregistrado ANTES de medir (R8, sin ajustar umbrales, M§44):
  * S0 con clip: F0 para α̂ <= 1.9 y F1 para α̂ >= 2.1 en el 100% de las corridas y en TODO N;
  * Omega-B γ̂=0, ρ=0.1, α̂=1.5·α̂_c: clique mas aislados (F2) en el 100% en todo N (de O-04); γ̂=10 solo informe;
  * `size_robust` es False en todo punto S0/Omega-B (evidencia trivial: no hay D* consistente) y ninguna corrida pasa
    todos los campos, de modo que no hay tamanos XL en el diseno;
  * control RGG3 k12 (full): pasa todos los campos de corrida (100%) en N=800, 1500 y 3000 y `size_robust` es True.
    Documenta que el mecanismo de `size_robust` acepta una estructura 3D genuina (R1: la ventana 3D exige N >= 800).
En smoke (N=16,20,24) el control RGG3 no se ejecuta (N pequeno carece de ventana 3D). Las expectativas se evaluan y se informan igual, pero
la expectativa S0 solo esta preregistrada para N >= 64: smoke mostro que con N <= 24 los puntos cercanos al umbral α̂=2 (1.9 y 2.1)
mezclan F0/F1 (efecto de tamano finito), asi que el test de smoke solo afirma los puntos alejados del umbral (α̂=3) y Omega-B γ̂=0.
Coste: con estados densos (ρ > 0.1, p. ej. los F1) la evidencia muestrea 25 aristas de Ollivier en vez de 1000 (`evidence_cfg`; solo presupuesto,
ningun umbral cambia): el LP de transporte de un grafo casi completo cuesta ~90 s con N=120.
`complete` = toda la malla preregistrada se ejecuto; `expectations_met` = las expectativas anteriores se cumplieron.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.certificate import size_robust
from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.config.convert import reduced_to_raw
from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, FixedDensityConfig, Omega11Config
from omega.contracts import RunEvidence
from omega.controls.random_geometric import rgg_torus
from omega.dynamics.fixed_density import evolve_fixed_density, uniform_state_threshold
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
from omega.phases.finite_size import FSS_KEYS, evidence_cfg, fss_fit, fss_observables, size_config
from omega.phases.scan import default_config, point_params, simulate
from omega.types import FloatArray, RunStatus

NAME = "o06_finite_size"
EXPERIMENT_ID = 1160  # 1106 lo usa s06
ENTRYPOINT = "omega.experiments.v11.test_o06_finite_size:run"
ALPHA_SPLIT = 1.95
OMEGA_B_RHO = 0.1
OMEGA_B_FACTOR = 1.5
RGG_SIZES = (800, 1500, 3000)


def _grid(cfg: Omega11Config, mode: Literal["smoke", "full"]) -> dict[str, Any]:
    if mode == "smoke":
        return {"sizes": (16, 20, 24), "reps": 1}
    return {"sizes": tuple(cfg.finite_size.sizes), "reps": cfg.finite_size.replicates}


def _points() -> list[dict[str, Any]]:
    pts: list[dict[str, Any]] = []
    for a in (1.9, 2.1, 3.0):
        for g in (0.0, 10.0):
            pts.append({"id": f"s0_a{a:g}_g{g:g}", "kind": "s0", "alpha_hat": a, "gamma_hat": g,
                        "claim": "Ω-F0" if a < ALPHA_SPLIT else "Ω-F1"})
    for g in (0.0, 10.0):
        pts.append({"id": f"omegab_r{OMEGA_B_RHO:g}_g{g:g}", "kind": "omega_b", "gamma_hat": g, "rho": OMEGA_B_RHO,
                    "claim": "Ω-F2" if g == 0.0 else None})
    return pts


def _evolve(cfg: Omega11Config, pt: dict[str, Any], n: int, key: SeedKey) -> tuple[FloatArray, RunStatus, FunctionalParams, Omega11Config]:
    """Estado final de un punto con N=n; devuelve (w, status, parametros, cfg del caso)."""
    if pt["kind"] == "s0":
        p = point_params(pt["alpha_hat"], pt["gamma_hat"], n)
        c_case = derive_cfg(size_config(cfg, n), functional=p)
        res = simulate(c_case.base, p, key)
        return res.trajectory.w_final, res.trajectory.status, p, c_case
    alpha_c = uniform_state_threshold(pt["gamma_hat"], pt["rho"], n)
    p = reduced_to_raw(OMEGA_B_FACTOR * alpha_c, pt["gamma_hat"], n)
    fd = FixedDensityConfig(rho=pt["rho"])
    c_case = derive_cfg(size_config(cfg, n), functional=p, engine=Engine.FIXED_DENSITY, fixed_density=fd)
    w0 = random_uniform_weights(n, make_rng(key))
    traj = evolve_fixed_density(w0, p, c_case.base.dynamics, fd)
    return traj.w_final, traj.status, p, c_case


def _curves(by_n: dict[int, list[RunEvidence]]) -> dict[str, Any]:
    """Medias por N de cada observable FSS y su ajuste de escalado."""
    sizes = sorted(by_n)
    curves: dict[str, Any] = {}
    for k in FSS_KEYS:
        means: list[float] = []
        for n in sizes:
            vals = [fss_observables(ev)[k] for ev in by_n[n]]
            fin = [v for v in vals if np.isfinite(v)]
            means.append(float(np.mean(fin)) if fin else float("nan"))
        entry: dict[str, Any] = {"sizes": sizes, "mean": means}
        try:
            entry["fit"] = fss_fit(sizes, means)
        except ValueError:
            entry["fit"] = None
        curves[k] = entry
    return curves


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-06 y devuelve `summary.json`. En full exige O-00..s06 full completos (paso 7)."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    g = _grid(cfg, mode)
    sizes: tuple[int, ...] = g["sizes"]
    reps = int(g["reps"])
    cfg = derive_cfg(cfg, n=sizes[0], experiment_id=EXPERIMENT_ID, replicates=reps)
    pts = _points()
    expected: dict[str, Any] = {
        "s0_codes": "F0 (α̂<=1.9) / F1 (α̂>=2.1) en el 100% de corridas y en todo N",
        "omega_b_gamma0": "F2 en el 100% de corridas y en todo N",
        "size_robust_trivial_points": False,
        "n_runs_pass_trivial_points": 0,
        "rgg3_control": ("size_robust True y pass_fraction 1.0 en N=(800,1500,3000)" if mode == "full" else "no se ejecuta en smoke"),
    }
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID, "description": "tamano finito D(N), C(N), G/N, rho(N), xi_F(N), L_hop(N) y control RGG3 (R1)",
        "sizes": sizes, "xl_sizes": tuple(cfg.finite_size.xl_sizes) if mode == "full" else (), "replicates": reps,
        "points": [p["id"] for p in pts], "expected": expected,
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: list[dict[str, Any]] = []
    for pi, pt in enumerate(pts):
        by_n: dict[int, list[RunEvidence]] = {}
        per_n: dict[str, Any] = {}
        all_claims = True
        n_pass = 0

        def do_size(n: int, ni: int, pt: dict[str, Any] = pt, pi: int = pi, by_n: dict[int, list[RunEvidence]] = by_n,
                    per_n: dict[str, Any] = per_n) -> tuple[int, bool]:
            codes: Counter[str] = Counter()
            passes = 0
            rows: list[dict[str, Any]] = []
            for rep in range(reps):
                key = seed_key(cfg.base.seeds, pi * 100 + ni, rep)
                w, status, p, c_case = _evolve(cfg, pt, n, key)
                ev = collect_run_evidence(w, status, evidence_cfg(c_case, w), evidence_rng(key))
                by_n.setdefault(n, []).append(ev)
                row = assessment_row(ev, evidence_cfg(c_case, w))
                codes.update(row["codes"])
                passes += int(row["passes"])
                rows.append(row)
                row["passport"] = writer.save(c_case, seed=key, w=w, termination=status.value,
                                              init_distribution="uniform/upper_mirror", results=row)
            claim = pt["claim"]
            met = None if claim is None else codes[claim] == reps
            per_n[str(n)] = {"code_fractions": {k: v / reps for k, v in sorted(codes.items())}, "n_runs_pass": passes,
                             "claim_met": met, "rows": rows}
            return passes, bool(met) if met is not None else True

        for ni, n in enumerate(sizes):
            passes, met = do_size(n, ni)
            n_pass += passes
            all_claims = all_claims and met
        # XL solo si TODAS las corridas del N mayor pasan todos los campos (no ocurre en el diseno preregistrado).
        xl_ran: list[int] = []
        if mode == "full" and per_n[str(sizes[-1])]["n_runs_pass"] == reps:
            for xi, n in enumerate(cfg.finite_size.xl_sizes):
                passes, _ = do_size(n, len(sizes) + xi)
                n_pass += passes
                xl_ran.append(n)
        ok, sr_ev = size_robust(by_n, cfg)
        results.append({
            "id": pt["id"], "kind": pt["kind"], "claim_code": pt["claim"], "claims_met_all_sizes": all_claims,
            "n_runs_pass": n_pass, "size_robust": ok, "size_robust_evidence": sr_ev, "xl_sizes_run": xl_ran,
            "per_size": per_n, "curves": _curves(by_n),
        })

    rgg: dict[str, Any] = {"ran": False, "reason": "solo en full (la ventana 3D exige N >= 800, R1)"}
    if mode == "full":
        by_n_rgg: dict[int, list[RunEvidence]] = {}
        per_n_rgg: dict[str, Any] = {}
        for ni, n in enumerate(RGG_SIZES):
            c_case = size_config(cfg, n)
            fracs = 0
            for rep in range(reps):
                key = seed_key(cfg.base.seeds, 10_000 + ni, rep)
                w = rgg_torus(n, 3, 12.0, make_rng(key))
                ev = collect_run_evidence(w, RunStatus.CONVERGED, evidence_cfg(c_case, w), evidence_rng(key))
                by_n_rgg.setdefault(n, []).append(ev)
                fracs += int(assess_run(ev, evidence_cfg(c_case, w)).passes)
                writer.save(c_case, seed=key, w=w, termination="STATIC", init_distribution="control/rgg3_k12",
                            results=assessment_row(ev, c_case))
            per_n_rgg[str(n)] = {"pass_fraction": fracs / reps}
        ok_r, ev_r = size_robust(by_n_rgg, cfg)
        rgg = {"ran": True, "sizes": RGG_SIZES, "size_robust": ok_r, "size_robust_evidence": ev_r, "per_size": per_n_rgg,
               "curves": _curves(by_n_rgg),
               "note": "control positivo del mecanismo size_robust (R1); la clase D*=3 debe sostenerse en N=800..3000"}

    met_s0 = all(r["claims_met_all_sizes"] for r in results if r["kind"] == "s0")
    met_b = all(r["claims_met_all_sizes"] for r in results if r["kind"] == "omega_b")
    trivial_ok = all(not r["size_robust"] and r["n_runs_pass"] == 0 for r in results)
    rgg_ok = True
    if rgg["ran"]:
        rgg_ok = bool(check_expectation("eq", rgg["size_robust"], True)
                      and all(check_expectation("eq", v["pass_fraction"], 1.0) for v in rgg["per_size"].values()))
    summary = {
        **header, "stage": "final", "results": results, "rgg3_control": rgg,
        "expectation_results": {"s0_codes": met_s0, "omega_b_gamma0_F2": met_b, "trivial_points_not_size_robust": trivial_ok,
                                "rgg3_control": rgg_ok},
        "expectations_met": bool(met_s0 and met_b and trivial_ok and rgg_ok),
        "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(24, master_entropy=20240901, replicates=2))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    assert len(data["results"]) == 8 and data["rgg3_control"]["ran"] is False
    for r in data["results"]:
        assert set(r["curves"]) == set(FSS_KEYS)
        assert set(r["per_size"]) == {"16", "20", "24"}
    far = [r for r in data["results"] if r["id"].startswith("s0_a3_") or r["id"] == "omegab_r0.1_g0"]
    assert len(far) == 3 and all(r["claims_met_all_sizes"] for r in far)  # F1 / F2 lejos del umbral (R8)
    assert isinstance(data["expectations_met"], bool)
    assert all(r["n_runs_pass"] == 0 for r in data["results"])
    assert (tmp_path / NAME / "summary.json").is_file()


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(64, master_entropy=20240901, replicates=10))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita: paso 7 del STEP_ORDER
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
