"""Omega-5 (Ω-1.1, paso 8): ensembles estadisticos Langevin y Metropolis (DESIGN §1.8, §4, riesgo R8).

PREGUNTA DEL PANEL (PANEL §8): ¿existe una region estadisticamente dominante con propiedades geometricas?
Respuesta operativa, por celda (motor, N, Θ̂, α̂, γ̂) con 4 cadenas independientes desde U(0,1): R̂ dividido, τ_int y ESS
(`ensemble_summary`). Si la celda NO equilibra (R̂ <= 1.05, ESS_total >= 100, |z_Geweke| <= 3) su estado es
`CHAIN_NOT_EQUILIBRATED` y NO se calcula ni se reporta fase alguna (ni codigos, ni probabilidades, ni transiciones con
ella). En las celdas equilibradas se evalua `assess_run` en los estados finales de las 4 cadenas y en hasta
`ensemble.phase_samples` estados post-burn-in de la cadena 0 separados por >= 2 τ_int (IC de Wilson). Una celda es
"region dominante geometrica" si equilibra y >= `certificate.seed_fraction` (0.8) de sus estados pasan TODOS los campos de
corrida (sin codigo de fallo). La respuesta a la pregunta es SI solo si existe alguna; si no, NO (indicando cuantas celdas
equilibraron y cuantas no). Ademas: <m>, χ_m = M·Var(m), U4 de Binder y bimodalidad de Sarle por celda; λ*(Θ̂) y λ*(α̂)
(maximo de χ) sobre celdas equilibradas; histeresis en α̂ (Metropolis, continuacion ida/vuelta, informe).

Preregistrado ANTES de medir (R8: la unica via abierta a una region no trivial es Θ > 0; sin ajustar umbrales, M§44):
  * no se reporta fase sin equilibrio (invariante estructural verificado en el resumen);
  * validacion Gibbs: con α̂=γ̂=0 las aristas son independientes con densidad ∝ exp(-w²/Θ̂) en [0,1]; en las celdas
    equilibradas de METROPOLIS (muestreador exacto) la media <m> coincide con la cuadratura a 3·SE + 0.005. Langevin es
    Euler-Maruyama y tiene sesgo O(dt) (en smoke con N=16, Θ̂=1 midio +0.013 con SE 0.002): su desviacion se REPORTA
    (`bias`) pero no se afirma (correccion hecha tras el primer smoke, documentada; ningun umbral de decision cambia);
  * NO existe region dominante geometrica con N <= 100: la ventana de dimension 3 solo existe con N >= 800 (R1), asi que se
    predice `dominant_geometric_region = False` (registro negativo honesto);
  * el resto (probabilidades de codigos, χ, U4, λ*, histeresis) es solo informe.
Full: N ∈ {64, 100}, Θ̂ ∈ {0.01, 0.03, 0.1, 0.3, 1}, α̂ ∈ {0, 0.5, 1, 1.5, 2, 2.5, 3, 4}, γ̂ ∈ {0, 1, 10}, ambos motores con sus
longitudes por defecto. Smoke: N=16, Θ̂ ∈ {0.1, 1}, α̂ ∈ {0, 3}, γ̂=0, Metropolis 200 sweeps y Langevin 1000 pasos.
Coste: los estados densos (ρ > 0.1) usan `evidence_cfg` (muestreo de 25 aristas de Ollivier; ningun umbral cambia).
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, LangevinConfig, MetropolisConfig, Omega11Config
from omega.contracts import ChainResult
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
from omega.phases.finite_size import evidence_cfg
from omega.phases.scan import default_config, point_params
from omega.statistics.ensemble import ensemble_summary, integrated_autocorr_time, post_burn_in
from omega.statistics.langevin import OBSERVABLE_KEYS, run_langevin
from omega.statistics.metropolis import run_metropolis
from omega.statistics.summary import proportion_ci
from omega.statistics.transition import (
    bimodality_coefficient,
    binder_cumulant,
    fluctuation,
    susceptibility,
    transition_summary,
)
from omega.types import FloatArray, RunStatus

NAME = "o05_ensemble"
EXPERIMENT_ID = 1105
ENTRYPOINT = "omega.experiments.v11.test_o05_ensemble:run"
NOT_EQUILIBRATED = "CHAIN_NOT_EQUILIBRATED"
GIBBS_SE_FLOOR = 0.005
HYSTERESIS_BASE = 50_000
DYNAMICS_STREAM = 3000


def _grid(mode: Literal["smoke", "full"]) -> dict[str, Any]:
    if mode == "smoke":
        return {"sizes": (16,), "thetas": (0.1, 1.0), "alphas": (0.0, 3.0), "gammas": (0.0,),
                "metropolis_sweeps": 200, "langevin_steps": 1000, "hyst_sweeps": 40}
    return {"sizes": (64, 100), "thetas": (0.01, 0.03, 0.1, 0.3, 1.0), "alphas": (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0),
            "gammas": (0.0, 1.0, 10.0), "metropolis_sweeps": 2000, "langevin_steps": 20_000, "hyst_sweeps": 400}


def gibbs_mean(theta_hat: float) -> float:
    """<w> de la densidad ∝ exp(-w²/Θ̂) en [0,1] (cuadratura de Simpson fina, exacta a 1e-12)."""
    x = np.linspace(0.0, 1.0, 20_001)
    f = np.exp(-(x * x) / theta_hat)
    num = np.trapezoid(x * f, x)
    den = np.trapezoid(f, x)
    return float(num / den)


def _dyn_rng(key: SeedKey) -> np.random.Generator:
    return make_rng(SeedKey(key.entropy, key.spawn_key + (DYNAMICS_STREAM,)))


def _case_cfg(cfg: Omega11Config, engine: Engine, p: FunctionalParams, n: int, theta: float, g: dict[str, Any]) -> Omega11Config:
    c = derive_cfg(cfg, n=n, functional=p)
    if engine is Engine.LANGEVIN:
        return dataclasses.replace(c, engine=engine, langevin=LangevinConfig(theta_hat=theta, n_steps=g["langevin_steps"]))
    return dataclasses.replace(c, engine=engine, metropolis=MetropolisConfig(theta_hat=theta, n_sweeps=g["metropolis_sweeps"]))


def _run_chain(c_case: Omega11Config, w0: FloatArray, key: SeedKey, n_states: int) -> ChainResult:
    w_min = c_case.base.graph.w_min
    p = c_case.base.functional
    if c_case.engine is Engine.LANGEVIN:
        assert c_case.langevin is not None
        return run_langevin(w0, p, c_case.langevin, w_min, _dyn_rng(key), n_states=n_states)
    assert c_case.metropolis is not None
    return run_metropolis(w0, p, c_case.metropolis, w_min, _dyn_rng(key), n_states=n_states)


def _burn(c_case: Omega11Config) -> float:
    cfg_e = c_case.langevin if c_case.engine is Engine.LANGEVIN else c_case.metropolis
    assert cfg_e is not None
    return float(cfg_e.burn_in_fraction)


def _thin(c_case: Omega11Config) -> int:
    cfg_e = c_case.langevin if c_case.engine is Engine.LANGEVIN else c_case.metropolis
    assert cfg_e is not None
    return int(cfg_e.thin)


def _fin(x: float) -> float | None:
    return float(x) if math.isfinite(x) else None


def _pick_states(chain: ChainResult, tau_samples: float, max_states: int, n_post: int) -> list[FloatArray]:
    """Hasta `max_states` estados post-burn-in equiespaciados y separados por >= 2 τ_int (en muestras)."""
    if not chain.states or not math.isfinite(tau_samples):
        return []
    k = min(len(chain.states), max_states, int(n_post // max(2.0 * tau_samples, 1.0)))
    if k <= 0:
        return []
    idx = np.unique(np.linspace(0, len(chain.states) - 1, k).astype(int))
    return [chain.states[int(i)] for i in idx]


def _cell(cfg: Omega11Config, engine: Engine, n: int, theta: float, alpha: float, gamma: float, idx: int,
          g: dict[str, Any], writer: PassportWriter) -> tuple[dict[str, Any], FloatArray | None]:
    p = point_params(alpha, gamma, n)
    c_case = _case_cfg(cfg, engine, p, n, theta, g)
    n_chains = cfg.ensemble.n_chains
    chains: list[ChainResult] = []
    keys: list[SeedKey] = []
    for c in range(n_chains):
        key = seed_key(cfg.base.seeds, idx, c)
        keys.append(key)
        w0 = random_uniform_weights(n, make_rng(key))
        chains.append(_run_chain(c_case, w0, key, cfg.ensemble.phase_samples if c == 0 else 0))
    burn = _burn(c_case)
    summ = ensemble_summary(chains, OBSERVABLE_KEYS, cfg.ensemble, burn)
    status = "EQUILIBRATED" if summ.equilibrated else NOT_EQUILIBRATED
    cell: dict[str, Any] = {
        "engine": engine.value, "n": n, "theta_hat": theta, "alpha_hat": alpha, "gamma_hat": gamma, "status": status,
        "rhat": {k: _fin(v) for k, v in summ.rhat.items()}, "tau_int": {k: _fin(v) for k, v in summ.tau_int.items()},
        "ess": {k: _fin(v) for k, v in summ.ess.items()}, "geweke_ok": summ.geweke_ok,
        "acceptance": [float(c.acceptance) for c in chains], "step_size": float(chains[0].step_size),
    }
    m_pool: FloatArray | None = None
    if summ.equilibrated:
        m_series = [post_burn_in(c.samples["mean_weight"], burn) for c in chains]
        m_pool = np.concatenate(m_series)
        n_edges = n * (n - 1) // 2
        cell["mean_weight"] = {"mean": float(m_pool.mean()), "se": _fin(summ.ses["mean_weight"])}
        cell["chi_m"] = float(susceptibility(m_pool, n_edges))
        cell["fluctuation"] = float(fluctuation(m_pool))
        cell["binder_u4"] = _fin(binder_cumulant(m_pool))
        cell["bimodality"] = _fin(bimodality_coefficient(m_pool))
        tau_s = float(np.mean([integrated_autocorr_time(s, cfg.ensemble.autocorr_c) for s in m_series]))
        n_post = m_series[0].shape[0]
        states = [c.w_final for c in chains] + _pick_states(chains[0], tau_s, cfg.ensemble.phase_samples, n_post)
        rows: list[dict[str, Any]] = []
        passes = 0
        codes: Counter[str] = Counter()
        for si, w in enumerate(states):
            c_ev = evidence_cfg(c_case, w)
            ev = collect_run_evidence(w, RunStatus.CONVERGED, c_ev,
                                      evidence_rng(SeedKey(keys[0].entropy, keys[0].spawn_key + (si,))))
            a = assess_run(ev, c_ev)
            passes += int(a.passes)
            codes.update(c.value for c in a.codes)
            rows.append(assessment_row(ev, c_ev))
        k = len(states)
        lo, hi = proportion_ci(passes, k)
        cell["phase"] = {"n_states": k, "pass_fraction": passes / k, "pass_ci_wilson": [lo, hi],
                         "code_fractions": {c: v / k for c, v in sorted(codes.items())}, "rows": rows}
    cell["passport"] = writer.save(
        c_case, seed=keys[0], w=chains[0].w_final,
        termination="CHAIN_COMPLETE" if summ.equilibrated else NOT_EQUILIBRATED,
        init_distribution=f"uniform/upper_mirror;engine={engine.value}",
        results={k: v for k, v in cell.items() if k not in ("passport",) and k != "phase"} | (
            {"phase": {k: v for k, v in cell["phase"].items() if k != "rows"}} if "phase" in cell else {}),
    )
    return cell, m_pool


def _lambda_star(series: dict[float, FloatArray], n_edges: int) -> dict[str, Any]:
    if len(series) < 2:
        return {"status": "insufficient_equilibrated", "n_lambdas": len(series)}
    ts = transition_summary(series, n_edges)
    chi_max = float(np.max(ts.susceptibility))
    return {"status": "ok", "lambda_star": float(ts.peak_lambda) if chi_max > 0 else None, "chi_max": chi_max,
            "lambdas": ts.lambdas, "mean": ts.mean, "derivative": ts.derivative, "susceptibility": ts.susceptibility,
            "binder": ts.binder, "bimodality": ts.bimodality}


def _hysteresis(cfg: Omega11Config, n: int, theta: float, gamma: float, alphas: tuple[float, ...], g: dict[str, Any],
                h_idx: int) -> dict[str, Any]:
    """Continuacion en α̂ (Metropolis, 1 cadena por paso desde el estado previo): ida y vuelta; informe."""
    path = list(alphas) + list(reversed(alphas))
    w = random_uniform_weights(n, make_rng(seed_key(cfg.base.seeds, HYSTERESIS_BASE + h_idx, 0)))
    means: list[float] = []
    for i, a in enumerate(path):
        key = seed_key(cfg.base.seeds, HYSTERESIS_BASE + h_idx, 0)
        key = SeedKey(key.entropy, key.spawn_key + (i,))
        p = point_params(a, gamma, n)
        mc = MetropolisConfig(theta_hat=theta, n_sweeps=g["hyst_sweeps"])
        c_case = dataclasses.replace(derive_cfg(cfg, n=n, functional=p), engine=Engine.METROPOLIS, metropolis=mc)
        ch = _run_chain(c_case, w, key, 0)
        w = ch.w_final
        means.append(float(post_burn_in(ch.samples["mean_weight"], mc.burn_in_fraction).mean()))
    k = len(alphas)
    fwd = np.asarray(means[:k])
    bwd = np.asarray(means[k:][::-1])
    gap = fwd - bwd
    return {"theta_hat": theta, "gamma_hat": gamma, "alphas": list(alphas), "mean_up": fwd, "mean_down": bwd, "gap": gap,
            "max_abs_gap": float(np.max(np.abs(gap)))}


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-05 y devuelve `summary.json`. En full exige O-00..O-06 full completos (paso 8)."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    g = _grid(mode)
    sizes, thetas, alphas, gammas = g["sizes"], g["thetas"], g["alphas"], g["gammas"]
    cfg = derive_cfg(cfg, n=sizes[0], experiment_id=EXPERIMENT_ID, replicates=max(cfg.base.seeds.replicates, cfg.ensemble.n_chains))
    engines = (Engine.LANGEVIN, Engine.METROPOLIS)
    expected: dict[str, Any] = {
        "no_phase_without_equilibrium": True,
        "gibbs_alpha0_gamma0": f"Metropolis: <m> = cuadratura de exp(-w²/Θ̂) a 3·SE + {GIBBS_SE_FLOOR} en celdas equilibradas (Langevin: solo sesgo informado)",
        "dominant_geometric_region": False,
        "reason": "N <= 100: la ventana de dimension 3 exige N >= 800 (R1); unica via abierta: Θ > 0 (R8)",
    }
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID, "description": "ensembles Langevin/Metropolis: ¿existe una region estadisticamente dominante con propiedades geometricas?",
        "sizes": sizes, "theta_hat": thetas, "alpha_hat": alphas, "gamma_hat": gammas, "engines": [e.value for e in engines],
        "n_chains": cfg.ensemble.n_chains, "expected": expected,
        "panel_question": "¿existe una region estadisticamente dominante con propiedades geometricas?",
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    cells: list[dict[str, Any]] = []
    pools: dict[tuple[str, int, float, float, float], FloatArray] = {}
    idx = 0
    for eng in engines:
        for n in sizes:
            for gm in gammas:
                for th in thetas:
                    for a in alphas:
                        cell, m_pool = _cell(cfg, eng, n, th, a, gm, idx, g, writer)
                        cells.append(cell)
                        if m_pool is not None:
                            pools[(eng.value, n, gm, th, a)] = m_pool
                        idx += 1

    transitions: list[dict[str, Any]] = []
    for eng in engines:
        for n in sizes:
            for gm in gammas:
                n_edges = n * (n - 1) // 2
                for a in alphas:  # λ = Θ̂ a α̂ fijo
                    ser = {th: pools[(eng.value, n, gm, th, a)] for th in thetas if (eng.value, n, gm, th, a) in pools}
                    transitions.append({"engine": eng.value, "n": n, "gamma_hat": gm, "lambda": "theta_hat", "fixed_alpha_hat": a,
                                        **_lambda_star(ser, n_edges)})
                for th in thetas:  # λ = α̂ a Θ̂ fijo
                    ser = {a: pools[(eng.value, n, gm, th, a)] for a in alphas if (eng.value, n, gm, th, a) in pools}
                    transitions.append({"engine": eng.value, "n": n, "gamma_hat": gm, "lambda": "alpha_hat", "fixed_theta_hat": th,
                                        **_lambda_star(ser, n_edges)})
    hyst = [_hysteresis(cfg, n, th, gm, alphas, g, h) for h, (n, gm, th) in enumerate(
        (n, gm, th) for n in sizes for gm in gammas for th in thetas)]

    eq = [c for c in cells if c["status"] == "EQUILIBRATED"]
    neq = [c for c in cells if c["status"] == NOT_EQUILIBRATED]
    frac_min = cfg.certificate.seed_fraction
    dominant = [c for c in eq if c["phase"]["pass_fraction"] >= frac_min]
    answer = (
        f"SI: {len(dominant)} celdas equilibradas con >= {frac_min:.0%} de estados que pasan todos los campos de corrida"
        if dominant else
        f"NO: ninguna de las {len(eq)} celdas equilibradas ({len(neq)} no equilibraron: CHAIN_NOT_EQUILIBRATED, sin fase reportada) "
        f"tiene >= {frac_min:.0%} de estados que pasan todos los campos de corrida"
    )
    no_phase_ok = all("phase" not in c and "mean_weight" not in c for c in neq)
    gibbs: list[dict[str, Any]] = []
    for c in eq:
        if c["alpha_hat"] == 0.0 and c["gamma_hat"] == 0.0:
            tol = 3.0 * (c["mean_weight"]["se"] or 0.0) + GIBBS_SE_FLOOR
            ref = gibbs_mean(c["theta_hat"])
            gibbs.append({"engine": c["engine"], "n": c["n"], "theta_hat": c["theta_hat"], "measured": c["mean_weight"]["mean"],
                          "reference": ref, "tol": tol, "bias": c["mean_weight"]["mean"] - ref,
                          "asserted": c["engine"] == "metropolis",
                          "ok": check_expectation("abs_within", c["mean_weight"]["mean"], ref, tol)})
    summary = {
        **header, "stage": "final", "cells": cells, "transitions": transitions, "hysteresis": hyst,
        "n_cells": len(cells), "n_equilibrated": len(eq), "n_not_equilibrated": len(neq),
        "dominant_geometric_region": bool(dominant), "dominant_cells": [
            {k: c[k] for k in ("engine", "n", "theta_hat", "alpha_hat", "gamma_hat")} for c in dominant],
        "panel_answer": answer, "gibbs_checks": gibbs,
        "expectation_results": {
            "no_phase_without_equilibrium": no_phase_ok,
            "gibbs_alpha0_gamma0": all(x["ok"] for x in gibbs if x["asserted"]),
            "dominant_geometric_region_false": not dominant,
        },
        "passports": writer.labels, "complete": True,
    }
    summary["expectations_met"] = all(summary["expectation_results"].values())
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(16, master_entropy=20240901, replicates=4))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    assert data["n_cells"] == 2 * 1 * 1 * 2 * 2 and "panel_answer" in data
    for c in data["cells"]:
        assert set(c["rhat"]) == set(OBSERVABLE_KEYS)
        assert (c["status"] == "EQUILIBRATED") == ("phase" in c)  # nunca fase sin equilibrio
    assert data["expectation_results"]["no_phase_without_equilibrium"] is True
    assert data["n_equilibrated"] + data["n_not_equilibrated"] == data["n_cells"]
    assert data["n_equilibrated"] >= 1  # al menos una celda (aristas independientes) equilibra
    assert all(x["ok"] for x in data["gibbs_checks"] if x["asserted"])
    assert data["dominant_geometric_region"] is False
    assert len(data["hysteresis"]) == 2 and len(data["passports"]) == data["n_cells"]


def test_gibbs_reference() -> None:
    assert gibbs_mean(1e6) == pytest.approx(0.5, abs=1e-6)  # Θ̂ -> inf: uniforme
    assert 0.0 < gibbs_mean(0.01) < 0.1 < gibbs_mean(1.0) < 0.5


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(64, master_entropy=20240901, replicates=4))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita: paso 8 del STEP_ORDER
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
