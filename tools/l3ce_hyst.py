"""L-3c-E'' (Rev. 3 §14): histeresis. Igual que L-3c-E' pero partiendo de la bicapa. Constantes congeladas.

Partes: n32 | n64 | n128_0 .. n128_4 | final
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from l3ce_prime import B_GRID, EVERY, MEAS, N128_SPLIT, STEP, THERM, b_c, classifier, flag, git_state, rng_for  # noqa: E402

from omega.causality.mcmc import run_chain  # noqa: E402

ENT_CHAIN_HYST = 20261028
UP = Path("results") / "l3ce_prime" / "summary.json"
OUT = Path("results") / "l3ce_hyst"


def bilayer_v(n: int) -> np.ndarray[Any, Any]:
    """u = identidad; dos antichains apiladas: v = [m-1..0, n-1..m]."""
    m = n // 2
    return np.concatenate([np.arange(m)[::-1], np.arange(m, n)[::-1]]).astype(np.int64)


def run_part(n: int, ks: list[int], therm: int, meas: int, every: int) -> dict[str, Any]:
    delta2, ref2 = classifier(n)
    rows = []
    for k in ks:
        b = B_GRID[k]
        rng = rng_for(ENT_CHAIN_HYST, 1000 * n + k, 0)
        res = run_chain(bilayer_v(n), b / n, therm, meas, every, rng)
        flags = [flag(c, delta2, ref2) for c in res.samples]
        rows.append({"k": k, "b": b, "flag_fraction": float(np.mean(flags)),
                     "mean_S": float(np.mean(res.actions)), "acceptance": res.acceptance})
    return {"N": n, "delta2": delta2, "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--part", default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "smoke":
        r = run_part(16, [0, 4, 8], 20, 20, 10)
        print([(round(x["b"], 2), x["flag_fraction"]) for x in r["rows"]])
        return
    if dirty:
        sys.exit("modo full exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["n32", "n64"] + [f"n128_{i}" for i in range(len(N128_SPLIT))]
    if args.part == "final":
        parts = {p: json.loads((OUT / f"{p}.json").read_text()) for p in names}
        commits = sorted({v["commit"] for v in parts.values()})
        rows = {32: parts["n32"]["data"]["rows"], 64: parts["n64"]["data"]["rows"],
                128: [r for i in range(len(N128_SPLIT)) for r in parts[f"n128_{i}"]["data"]["rows"]]}
        down = {n: b_c(r) for n, r in rows.items()}
        floor = {n: min(r, key=lambda x: x["b"])["flag_fraction"] < 0.5 for n, r in rows.items()}
        up = {int(k): float(v) for k, v in json.loads(UP.read_text())["b_c"].items()}
        gap = {n: round(abs(math.log(up[n] / down[n]) / math.log(STEP))) for n in rows}
        if any(floor.values()) or any(math.isinf(x) for x in down.values()):
            verdict = "INDETERMINADO"
        elif all(g <= 1 for g in gap.values()):
            verdict = "SIN HISTERESIS RELEVANTE"
        elif gap[128] >= 2 and gap[128] > gap[32]:
            verdict = "HISTERESIS CRECIENTE"
        else:
            ok = down[64] >= down[32] / STEP - 1e-9 and down[128] >= down[64] / STEP - 1e-9
            verdict = "HISTERESIS ESTABLE (PASA se mantiene)" if ok else "HISTERESIS ESTABLE (PASA se rebaja)"
        summary = {"phase": "final", "commits_partes": commits, "b_c_up": up, "b_c_down": down,
                   "gap_steps": gap, "floor": floor, "L-3c-E''": verdict}
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
        print(json.dumps(summary, indent=2, default=str))
        return
    t0 = time.time()
    if args.part == "n32":
        data = run_part(32, list(range(9)), THERM, MEAS, EVERY)
    elif args.part == "n64":
        data = run_part(64, list(range(9)), THERM, MEAS, EVERY)
    elif args.part is not None and args.part.startswith("n128_"):
        data = run_part(128, N128_SPLIT[int(args.part.split("_")[1])], THERM, MEAS, EVERY)
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
