"""L-3c-E' (Rev. 3 §12): estabilidad euclidea con acoplamiento beta = b/N. Constantes congeladas. No ajustar.

Partes: n32 | n64 | n128_0 .. n128_4 | final
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import interval_abundances, mm_dimension  # noqa: E402
from omega.causality.mcmc import run_chain  # noqa: E402
from omega.causality.orders import layered_2order, random_2order  # noqa: E402
from omega.causality.profile import K_MAX, d_ref, profile_distance, reference_abundances  # noqa: E402

EXPERIMENT_ID = 304
ENT_CHAIN, ENT_CAL, ENT_VAL = 20261025, 20261026, 20261027
B_GRID = [2.0 * 2 ** (k / 4) for k in range(9)]
STEP = 2 ** (1 / 4)
THERM, MEAS, EVERY, N_CAL = 1000, 1000, 50, 40
N128_SPLIT = [[0, 1], [2, 3], [4, 5], [6, 7], [8]]
OUT = Path("results") / "l3ce_prime"


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def classifier(n: int) -> tuple[float, np.ndarray[Any, Any]]:
    ref = reference_abundances(2, n, K_MAX)
    d = [profile_distance(interval_abundances(random_2order(n, rng_for(ENT_CAL, n, s)), K_MAX), ref) for s in range(N_CAL)]
    return 1.5 * float(np.percentile(d, 95)), ref


def flag(c: np.ndarray[Any, Any], delta2: float, ref2: np.ndarray[Any, Any]) -> bool:
    dr = d_ref(mm_dimension(c))
    if dr is None or dr > 2:
        return False
    ref = ref2 if dr == 2 else reference_abundances(1, c.shape[0], K_MAX)
    return bool(profile_distance(interval_abundances(c, K_MAX), ref) <= delta2)


def validate(n: int, delta2: float, ref2: np.ndarray[Any, Any]) -> dict[str, Any]:
    pos = float(np.mean([flag(random_2order(n, rng_for(ENT_VAL, n, s)), delta2, ref2) for s in range(100)]))
    lay = {
        str(k): float(np.mean([flag(layered_2order(n, k, rng_for(ENT_VAL, 100 * n + k, s)), delta2, ref2) for s in range(20)]))
        for k in (2, 4, 8)
    }
    return {"pos": pos, "layered": lay, "ok": pos >= 0.9 and all(v <= 0.1 for v in lay.values())}


def run_part(n: int, ks: list[int], therm: int, meas: int, every: int, with_validation: bool) -> dict[str, Any]:
    delta2, ref2 = classifier(n)
    out: dict[str, Any] = {"N": n, "delta2": delta2, "rows": []}
    if with_validation:
        out["validation"] = validate(n, delta2, ref2)
    for k in ks:
        b = B_GRID[k]
        rng = rng_for(ENT_CHAIN, 1000 * n + k, 0)
        v0 = rng.permutation(n).astype(np.int64)
        res = run_chain(v0, b / n, therm, meas, every, rng)
        flags = [flag(c, delta2, ref2) for c in res.samples]
        out["rows"].append({"k": k, "b": b, "beta": b / n, "flag_fraction": float(np.mean(flags)),
                            "mean_S": float(np.mean(res.actions)), "acceptance": res.acceptance})
    return out


def b_c(rows: list[dict[str, Any]]) -> float:
    for r in sorted(rows, key=lambda x: x["b"]):
        if r["flag_fraction"] < 0.5:
            return float(r["b"])
    return float("inf")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--part", default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "smoke":
        t0 = time.time()
        r = run_part(16, [0, 4, 8], 20, 20, 10, True)
        print(json.dumps({"validation": r["validation"], "rows": [(round(x["b"], 2), x["flag_fraction"]) for x in r["rows"]],
                          "s": round(time.time() - t0, 1)}))
        return
    if dirty:
        sys.exit("modo full exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["n32", "n64"] + [f"n128_{i}" for i in range(len(N128_SPLIT))]
    if args.part == "final":
        parts = {p: json.loads((OUT / f"{p}.json").read_text()) for p in names}
        commits = sorted({v["commit"] for v in parts.values()})
        rows128 = [r for i in range(len(N128_SPLIT)) for r in parts[f"n128_{i}"]["data"]["rows"]]
        bc = {32: b_c(parts["n32"]["data"]["rows"]), 64: b_c(parts["n64"]["data"]["rows"]), 128: b_c(rows128)}
        val = {32: parts["n32"]["data"]["validation"], 64: parts["n64"]["data"]["validation"], 128: parts["n128_0"]["data"]["validation"]}
        if not all(v["ok"] for v in val.values()) or any(x == float("inf") for x in bc.values()):
            verdict = "INDETERMINADO"
        elif bc[64] >= bc[32] / STEP - 1e-9 and bc[128] >= bc[64] / STEP - 1e-9:
            verdict = "PASA (ventana estable con beta = b/N)"
        else:
            verdict = "FALLA"
        summary = {"phase": "final", "commits_partes": commits, "un_solo_commit": len(commits) == 1,
                   "b_c": {str(k): v for k, v in bc.items()}, "validacion": {str(k): v for k, v in val.items()}, "L-3c-E'": verdict}
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    t0 = time.time()
    if args.part == "n32":
        data = run_part(32, list(range(9)), THERM, MEAS, EVERY, True)
    elif args.part == "n64":
        data = run_part(64, list(range(9)), THERM, MEAS, EVERY, True)
    elif args.part is not None and args.part.startswith("n128_"):
        i = int(args.part.split("_")[1])
        data = run_part(128, N128_SPLIT[i], THERM, MEAS, EVERY, i == 0)
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
