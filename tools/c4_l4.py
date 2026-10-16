"""C4 y L4 (Rev. 3 §19): validacion de la accion 4D y compatibilidad de la rama lorentziana en 4D.

Partes: c4 | l4_small | l4_1024 | l4_2048 | final. Constantes congeladas; no ajustar.
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

from omega.causality.estimators import between_counts, smearing_weight_4d  # noqa: E402
from omega.causality.reference4d import expected_action_4d, pair_volume_fractions  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 308
ENT_REF, ENT_C4, ENT_L4 = 20261038, 20261039, 20261040
REF_PAIRS = 4_000_000
C4_N, C4_EPS, C4_M = (128, 512), (1.0, 0.25, 0.05), 200
L4_N, L4_M = (128, 256, 512, 1024, 2048), 100
L4_FIXED = (1.0, 0.5, 0.25, 0.1, 0.05)
L4_K = (16, 64)
D_MAX = 3 / (2 * math.pi)
OUT = Path("results") / "c4_l4"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def s4_from_between(vals: np.ndarray[Any, Any], n: int, eps: float) -> float:
    if eps == 1.0:
        k = np.bincount(vals.astype(np.int64), minlength=4)[:4]
        return float(n - k[0] + 9 * k[1] - 16 * k[2] + 8 * k[3])
    return float(n - eps * float(np.sum(smearing_weight_4d(vals.astype(np.float64), eps))))


def schedules(n: int) -> dict[str, float]:
    out = {f"fijo_{e}": e for e in L4_FIXED}
    out.update({f"K{k}": k / n for k in L4_K})
    return out


def part_c4() -> dict[str, Any]:
    v = pair_volume_fractions(REF_PAIRS, rng_for(ENT_REF, 0))
    rows = {}
    for i, n in enumerate(C4_N):
        samples: dict[float, list[float]] = {e: [] for e in C4_EPS}
        for s in range(C4_M):
            b = between_counts(causal_matrix(sprinkle_diamond(n, 4, rng_for(ENT_C4, i, s))))
            vals = b[b >= 0]
            for e in C4_EPS:
                samples[e].append(s4_from_between(vals, n, e))
        for e in C4_EPS:
            ref, ref_se = expected_action_4d(n, e, v)
            mean = float(np.mean(samples[e]))
            se = float(np.std(samples[e], ddof=1) / math.sqrt(C4_M))
            tol = 3 * math.sqrt(se**2 + ref_se**2)
            rows[f"N{n}_eps{e}"] = {"mean": mean, "se": se, "ref": ref, "ref_se": ref_se, "ok": abs(mean - ref) <= tol}
    return {"rows": rows, "pass": all(r["ok"] for r in rows.values())}


def part_l4(ns: tuple[int, ...]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for n in ns:
        i = L4_N.index(n)
        sched = schedules(n)
        acc: dict[str, list[float]] = {k: [] for k in sched}
        for s in range(L4_M):
            b = between_counts(causal_matrix(sprinkle_diamond(n, 4, rng_for(ENT_L4, i, s))))
            vals = b[b >= 0]
            for k, e in sched.items():
                acc[k].append(s4_from_between(vals, n, e))
        out[str(n)] = {k: {"eps": sched[k], "sd": float(np.std(v, ddof=1)), "D4": float(np.std(v, ddof=1)) / sched[k]} for k, v in acc.items()}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["c4", "l4_small", "l4_1024", "l4_2048", "final"], required=True)
    args = ap.parse_args()
    head, dirty = git_state()
    if dirty:
        sys.exit("exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "final":
        c4 = json.loads((OUT / "c4.json").read_text())["data"]
        l4: dict[str, Any] = {}
        for p in ("l4_small", "l4_1024", "l4_2048"):
            l4.update(json.loads((OUT / f"{p}.json").read_text())["data"])
        fits = {}
        for key in schedules(128):
            x = np.log([float(n) for n in L4_N])
            y = np.log([l4[str(n)][key]["D4"] for n in L4_N])
            a = np.vstack([x, np.ones_like(x)]).T
            coef = np.linalg.lstsq(a, y, rcond=None)[0]
            r = y - a @ coef
            se = math.sqrt(float(r @ r) / (len(x) - 2) * float(np.linalg.inv(a.T @ a)[0, 0]))
            fits[key] = {"eta": float(coef[0]), "se": se, "D4_2048": l4["2048"][key]["D4"]}
        compat = [k for k, f in fits.items() if f["D4_2048"] <= D_MAX and f["eta"] + 2 * f["se"] <= 0]
        incompat = all(f["eta"] - 2 * f["se"] > 0 for f in fits.values())
        if not c4["pass"]:
            verdict = "INDETERMINADO (C4 no aprobado)"
        else:
            verdict = "COMPATIBLE" if compat else ("INCOMPATIBLE (esquemas probados)" if incompat else "INDETERMINADO")
        summary = {"commit": head, "C4": c4["pass"], "fits": fits, "L4": verdict, "compatibles": compat}
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    if args.part == "c4":
        data: dict[str, Any] = part_c4()
    elif args.part == "l4_small":
        data = part_l4((128, 256, 512))
    elif args.part == "l4_1024":
        data = part_l4((1024,))
    else:
        data = part_l4((2048,))
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s", data.get("pass", ""))


if __name__ == "__main__":
    main()
