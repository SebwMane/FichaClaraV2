"""C1 y L-3c-E (Rev. 3 §11). Constantes congeladas en docs/OMEGA_2_PROPUESTA_CAUSAL.md §11. No ajustar.

Por limite de tiempo del entorno, la corrida full se ejecuta por partes que escriben JSON parciales:
  --part c1 | l3ce32 | l3ce64a | l3ce64b | final
'final' verifica que todas las partes se hicieron en el mismo commit limpio y calcula los veredictos.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import bd_action_2d_smeared, interval_abundances, mm_dimension  # noqa: E402
from omega.causality.exact import expected_smeared_action_2order  # noqa: E402
from omega.causality.mcmc import run_chain  # noqa: E402
from omega.causality.orders import layered_2order, random_2order  # noqa: E402
from omega.causality.profile import K_MAX, d_ref, profile_distance, reference_abundances  # noqa: E402

EXPERIMENT_ID = 303
ENT_C1, ENT_CHAIN, ENT_CAL = 20261021, 20261022, 20261023
ENT_VAL = 20261024  # AUDITORIA A-14 / enmienda L3cE-A1
C1_N, C1_EPS, C1_M = (64, 256), (0.5, 0.25, 0.1), 200
BETAS = [0.025 * math.sqrt(2) ** k for k in range(13)]
THERM, MEAS, EVERY = 1000, 1000, 50
N_CAL = 40
OUT = Path("results") / "c1_l3ce"


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def part_c1(scale: float = 1.0) -> dict[str, Any]:
    out: dict[str, Any] = {}
    m = max(2, int(C1_M * scale))
    for i, n in enumerate(C1_N):
        orders = [random_2order(n, rng_for(ENT_C1, i, s)) for s in range(m)]
        for eps in C1_EPS:
            vals = [bd_action_2d_smeared(c, eps) for c in orders]
            mean, se = float(np.mean(vals)), float(np.std(vals, ddof=1) / math.sqrt(m))
            exact = expected_smeared_action_2order(n, eps)
            out[f"N{n}_eps{eps}"] = {"mean": mean, "se": se, "exact": exact, "ok": abs(mean - exact) <= 3 * se}
    return out


def calibrate_delta(n: int) -> tuple[float, np.ndarray[Any, Any]]:
    ref = reference_abundances(2, n, K_MAX)
    d = [profile_distance(interval_abundances(random_2order(n, rng_for(ENT_CAL, n, s)), K_MAX), ref) for s in range(N_CAL)]
    return 1.5 * float(np.percentile(d, 95)), ref


def flag(c: np.ndarray[Any, Any], delta2: float, ref2: np.ndarray[Any, Any]) -> bool:
    dr = d_ref(mm_dimension(c))
    if dr is None or dr > 2:
        return False  # en L-3c-E solo se calibra d = 2; d_ref 3 o 4 no aparece en 2-ordenes
    ref = ref2 if dr == 2 else reference_abundances(1, c.shape[0], K_MAX)
    return bool(profile_distance(interval_abundances(c, K_MAX), ref) <= delta2)


def validate_classifier(n: int, delta2: float, ref2: np.ndarray[Any, Any]) -> dict[str, Any]:
    """Enmienda L3cE-A1: el clasificador recalibrado debe separar 2-ordenes aleatorios de 2-ordenes en capas a este N."""
    pos = float(np.mean([flag(random_2order(n, rng_for(ENT_VAL, n, s)), delta2, ref2) for s in range(100)]))
    lay = {
        str(k): float(np.mean([flag(layered_2order(n, k, rng_for(ENT_VAL, 100 * n + k, s)), delta2, ref2) for s in range(20)]))
        for k in (2, 4, 8)
    }
    return {"pos": pos, "layered": lay, "ok": pos >= 0.9 and all(v <= 0.1 for v in lay.values())}


def part_l3ce(n: int, betas: list[float], therm: int, meas: int, every: int, k0: int = 0) -> dict[str, Any]:
    delta2, ref2 = calibrate_delta(n)
    validation = validate_classifier(n, delta2, ref2)
    rows = []
    for k, beta in enumerate(betas, start=k0):
        rng = rng_for(ENT_CHAIN, 1000 * n + k, 0)
        v0 = rng.permutation(n).astype(np.int64)
        res = run_chain(v0, beta, therm, meas, every, rng)
        flags = [flag(c, delta2, ref2) for c in res.samples]
        rows.append({
            "beta": beta, "flag_fraction": float(np.mean(flags)), "mean_S": float(np.mean(res.actions)),
            "acceptance": res.acceptance, "actions": res.actions,
        })
    return {"N": n, "delta2": delta2, "validation": validation, "rows": rows}


def beta_c(rows: list[dict[str, Any]]) -> float:
    for r in sorted(rows, key=lambda x: x["beta"]):
        if r["flag_fraction"] < 0.5:
            return float(r["beta"])
    return float("inf")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--part", choices=["c1", "l3ce32", "l3ce64a", "l3ce64b", "final"], default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "smoke":
        t0 = time.time()
        c1 = part_c1(scale=0.1)
        l3 = part_l3ce(16, BETAS[::4], 20, 20, 10)
        print(json.dumps({"c1": c1, "l3ce16": {k: v for k, v in l3.items() if k != "rows"},
                          "flags": [(round(r["beta"], 3), r["flag_fraction"]) for r in l3["rows"]],
                          "s": round(time.time() - t0, 1)}, indent=1))
        return
    if dirty:
        sys.exit("modo full exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.part == "final":
        parts = {p: json.loads((OUT / f"{p}.json").read_text()) for p in ("c1", "l3ce32", "l3ce64a", "l3ce64b")}
        commits = {v["commit"] for v in parts.values()}
        l64 = {"N": 64, "rows": parts["l3ce64a"]["data"]["rows"] + parts["l3ce64b"]["data"]["rows"]}
        bc32, bc64 = beta_c(parts["l3ce32"]["data"]["rows"]), beta_c(l64["rows"])
        summary = {
            "phase": "final", "commits_partes": sorted(commits), "un_solo_commit": len(commits) == 1,
            "C1-2": all(v["ok"] for v in parts["c1"]["data"].values()),
            "beta_c_32": bc32, "beta_c_64": bc64,
            "validacion_clasificador": {p: parts[p]["data"]["validation"] for p in ("l3ce32", "l3ce64a")},
            "L-3c-E": (
                "INDETERMINADO (clasificador no validado)"
                if not (parts["l3ce32"]["data"]["validation"]["ok"] and parts["l3ce64a"]["data"]["validation"]["ok"])
                else ("ventana estable" if bc64 >= bc32 else "ventana que se angosta")
            ),
        }
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    (OUT / "summary.json").write_text(json.dumps({"phase": "preregistered", "commit": head}, indent=2))
    t0 = time.time()
    if args.part == "c1":
        data: dict[str, Any] = part_c1()
    elif args.part == "l3ce32":
        data = part_l3ce(32, BETAS, THERM, MEAS, EVERY)
    elif args.part == "l3ce64a":
        data = part_l3ce(64, BETAS[:7], THERM, MEAS, EVERY)
    else:
        data = part_l3ce(64, BETAS[7:], THERM, MEAS, EVERY, k0=7)
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
