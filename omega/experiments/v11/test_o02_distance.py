"""Omega-2 (Ω-1.1): sensibilidad a la definicion de distancia por familia de estados (DESIGN §1.1, §4).

Se informa, por familia, la fraccion `distance_sensitive` y la dispersion `metric_spread` (HOP/INV/LOG). Preregistrado:
  * estados finales binarios de Omega-1.0 (S0): dispersion 0 por construccion cuando es finita (informe);
  * toros con pesos aleatorios: ninguno sensible (fraccion 0);
  * RGG euclidiano: sensible (fraccion 1; sensibilidad real de esa ponderacion con N=800, documentada);
  * el resto (sigmoide con max_steps acotado, U(0,1) estatico con w_min de sensibilidad, Langevin): solo informe, sin umbral.
`metric_spread` infinito (algun modo sin ventana) cuenta como sensible, como define §1.1.
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings11 import Engine, LangevinConfig, Omega11Config
from omega.controls.random_geometric import rgg_torus
from omega.experiments.reference_graphs import periodic_lattice
from omega.experiments.v11.gate import PassportWriter, derive_cfg, require_prerequisites, runs_root, write_summary
from omega.geometry.distance_suite import distance_suite
from omega.network.initialization import random_uniform_weights
from omega.phases.scan import default_config, point_params, simulate
from omega.statistics.langevin import run_langevin
from omega.types import FloatArray

NAME = "o02_distance"
EXPERIMENT_ID = 1102
ENTRYPOINT = "omega.experiments.v11.test_o02_distance:run"

# (state, termination, w0 | None, cfg_case)
Maker = Callable[[SeedKey, Omega11Config], tuple[FloatArray, str, FloatArray | None, Omega11Config]]


def _sym(w: FloatArray) -> FloatArray:
    u = np.triu(w, 1)
    return np.asarray(u + u.T, dtype=np.float64)


def _families(cfg: Omega11Config, mode: Literal["smoke", "full"]) -> list[dict[str, Any]]:
    n = 40 if mode == "smoke" else 200
    n_big = 40 if mode == "smoke" else 800
    steps = 200 if mode == "smoke" else 2000
    reps = 2 if mode == "smoke" else 5
    sig_steps = 1500 if mode == "smoke" else 20_000
    p_s0 = point_params(1.5, 1.0, n)
    fams: list[dict[str, Any]] = []

    def s0(sigmoid: bool) -> Maker:
        def make(key: SeedKey, c: Omega11Config) -> tuple[FloatArray, str, FloatArray | None, Omega11Config]:
            base = c.base
            if sigmoid:
                # El integrador sigmoide converge muy lentamente (MAX_STEPS en la practica): se acota a `sig_steps` pasos.
                base = dataclasses.replace(
                    base, dynamics=dataclasses.replace(base.dynamics, integrator="sigmoid", max_steps=sig_steps))
            res = simulate(base, p_s0, key)
            return res.trajectory.w_final, res.trajectory.status.value, res.w0, derive_cfg(c, functional=p_s0)
        return make

    def static(wm: float) -> Maker:
        def make(key: SeedKey, c: Omega11Config) -> tuple[FloatArray, str, FloatArray | None, Omega11Config]:
            return random_uniform_weights(c.base.init.n, make_rng(key)), "STATIC", None, derive_cfg(c, w_min=wm)
        return make

    def langevin(theta_hat: float) -> Maker:
        def make(key: SeedKey, c: Omega11Config) -> tuple[FloatArray, str, FloatArray | None, Omega11Config]:
            rng = make_rng(key)
            w0 = random_uniform_weights(c.base.init.n, rng)
            lc = LangevinConfig(theta_hat=theta_hat, n_steps=steps, thin=max(1, steps // 20))
            chain = run_langevin(w0, p_s0, lc, c.base.graph.w_min, rng)
            cc = dataclasses.replace(derive_cfg(c, functional=p_s0), engine=Engine.LANGEVIN, langevin=lc)
            return chain.w_final, "CHAIN_COMPLETE", w0, cc
        return make

    def torus_weighted(shape: tuple[int, ...]) -> Maker:
        def make(key: SeedKey, c: Omega11Config) -> tuple[FloatArray, str, FloatArray | None, Omega11Config]:
            rng = make_rng(key)
            lat = periodic_lattice(shape)
            w = lat * _sym(0.5 + 0.5 * rng.random(lat.shape))
            return np.asarray(w, dtype=np.float64), "STATIC", None, derive_cfg(c, n=lat.shape[0])
        return make

    def rgg_euclid(key: SeedKey, c: Omega11Config) -> tuple[FloatArray, str, FloatArray | None, Omega11Config]:
        w = rgg_torus(n_big, 3, 12.0, make_rng(key), "euclidean")
        return w, "STATIC", None, derive_cfg(c, n=n_big)

    fams += [
        {"name": "omega10_final_binary", "make": s0(False), "n": n, "reps": reps, "claim": ("max_finite_spread", 0.0)},
        {"name": "sigmoid_final", "make": s0(True), "n": n, "reps": reps, "claim": None},
    ]
    for wm in cfg.base.graph.w_min_sensitivity:
        fams.append({"name": f"static_uniform_wmin{wm:g}", "make": static(wm), "n": n, "reps": reps, "claim": None})
    for th in (0.1, 0.3):
        fams.append({"name": f"langevin_theta{th:g}", "make": langevin(th), "n": n, "reps": min(reps, 3), "claim": None})
    shapes = [(6, 6)] if mode == "smoke" else [(9, 9, 9), (28, 28)]
    for shp in shapes:
        fams.append({"name": "torus_random_weights_" + "x".join(map(str, shp)), "make": torus_weighted(shp),
                     "n": int(np.prod(shp)), "reps": min(reps, 3), "claim": ("fraction_sensitive", 0.0)})
    fams.append({"name": "rgg_euclidean_k12", "make": rgg_euclid, "n": n_big, "reps": min(reps, 3),
                 "claim": ("fraction_sensitive", 1.0)})
    return fams


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-02 y devuelve `summary.json`. En full exige O-00 y O-01 full completos."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    fams = _families(cfg, mode)
    max_reps = max(f["reps"] for f in fams)
    cfg = derive_cfg(cfg, experiment_id=EXPERIMENT_ID, replicates=max_reps)
    expected = [{"family": f["name"], "metric": f["claim"][0], "target": f["claim"][1]} for f in fams if f["claim"]]
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID,
        "description": "distance_sensitive y metric_spread por familia (DESIGN §1.1)",
        "expected": expected,
        "info_only_families": [f["name"] for f in fams if not f["claim"]],
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    for fi, fam in enumerate(fams):
        rows: list[dict[str, Any]] = []
        for rep in range(int(fam["reps"])):
            key = seed_key(cfg.base.seeds, fi, rep)
            fcfg = derive_cfg(cfg, n=int(fam["n"]))
            w, term, w0, c_case = fam["make"](key, fcfg)
            suite = distance_suite(w, c_case.base.graph, c_case.base.dimension, c_case.distance)
            row = {
                "distance_sensitive": bool(suite.distance_sensitive),
                "metric_spread": suite.metric_spread if math.isfinite(suite.metric_spread) else None,
                "resistance_exponent": suite.resistance_exponent if math.isfinite(suite.resistance_exponent) else None,
                "d": {m.value: (e.value if math.isfinite(e.value) else None) for m, e in suite.estimates.items()},
                "status": {m.value: e.status for m, e in suite.estimates.items()},
                "reason": suite.reason,
            }
            if rep == 0:
                row["passport"] = writer.save(c_case, seed=key, w=w, termination=term,
                                              init_distribution=f"family/{fam['name']}", results=row, w0=w0)
            rows.append(row)
        finite = [r["metric_spread"] for r in rows if r["metric_spread"] is not None]
        res: dict[str, Any] = {
            "n": int(fam["n"]), "n_states": len(rows),
            "fraction_sensitive": sum(1 for r in rows if r["distance_sensitive"]) / len(rows),
            "n_finite_spread": len(finite),
            "mean_finite_spread": float(np.mean(finite)) if finite else None,
            "max_finite_spread": float(np.max(finite)) if finite else None,
            "rows": rows,
        }
        if fam["claim"]:
            metric, target = fam["claim"]
            val = res[metric]
            # "max_finite_spread == 0" se cumple tambien si ninguna dispersion es finita (no hay ventana que comparar).
            res["claim_met"] = bool(val is None or val == target) if metric == "max_finite_spread" else bool(val == target)
        results[fam["name"]] = res
    summary = {
        **header, "stage": "final", "results": results,
        "expectations_met": all(r["claim_met"] for r in results.values() if "claim_met" in r),
        "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(40, master_entropy=20240901, replicates=2))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert "omega10_final_binary" in res and "rgg_euclidean_k12" in res and "langevin_theta0.1" in res
    assert all(0.0 <= r["fraction_sensitive"] <= 1.0 and r["n_states"] >= 2 for r in res.values())
    assert len(data["passports"]) == len(res)
    assert res["omega10_final_binary"]["claim_met"] is True  # binarios: dispersion 0 o sin ventana


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=5))
    data = json.loads(run(cfg, runs_root(), mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
