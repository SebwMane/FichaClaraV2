"""C4b y L4b (Rev. 3 §20). Partes: c4b | small | n2048 | final. Constantes congeladas."""

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

from omega.causality.reference4d import expected_action_4d  # noqa: E402
from omega.causality.torus import SLAB_T, retarded_action_4d, slab_causal_matrix, slab_pair_fractions, sprinkle_slab  # noqa: E402

EXPERIMENT_ID = 309
ENT_REF, ENT_C4B, ENT_L4B = 20261041, 20261042, 20261043
NS = (256, 512, 1024, 2048)
M = 100
FIXED = (1.0, 0.25, 0.05)
KS = (16, 64)
D_MAX = 3 / (2 * math.pi)
OUT = Path("results") / "l4b"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def schedules(n: int) -> dict[str, float]:
    out = {f"fijo_{e}": e for e in FIXED}
    out.update({f"K{k}": k / n for k in KS})
    return out


def part_c4b() -> dict[str, Any]:
    v = slab_pair_fractions(4_000_000, rng_for(ENT_REF, 0))
    n = 512
    rows = {}
    samples: dict[float, list[float]] = {0.25: [], 0.05: []}
    for s in range(200):
        c = slab_causal_matrix(sprinkle_slab(n, rng_for(ENT_C4B, s)))
        for e in samples:
            samples[e].append(float(retarded_action_4d(c, e).sum()))
    for e, vals in samples.items():
        ref, ref_se = expected_action_4d(n, e, v)
        mean, se = float(np.mean(vals)), float(np.std(vals, ddof=1) / math.sqrt(len(vals)))
        rows[str(e)] = {"mean": mean, "se": se, "ref": ref, "ref_se": ref_se, "ok": abs(mean - ref) <= 3 * math.sqrt(se**2 + ref_se**2)}
    return {"rows": rows, "pass": all(r["ok"] for r in rows.values())}


def part_l4b(ns: tuple[int, ...]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for n in ns:
        i = NS.index(n)
        sch = schedules(n)
        tot: dict[str, list[float]] = {k: [] for k in sch}
        intr: dict[str, list[float]] = {k: [] for k in sch}
        n_int: list[int] = []
        for s in range(M):
            x = sprinkle_slab(n, rng_for(ENT_L4B, i, s))
            c = slab_causal_matrix(x)
            top = x[:, 0] >= SLAB_T / 2
            n_int.append(int(top.sum()))
            for k, e in sch.items():
                r = retarded_action_4d(c, e)
                tot[k].append(float(r.sum()))
                intr[k].append(float(r[top].sum()))
        out[str(n)] = {k: {"eps": sch[k], "D_tot": float(np.std(tot[k], ddof=1)) / sch[k],
                           "D_int": float(np.std(intr[k], ddof=1)) / sch[k]} for k in sch}
        out[str(n)]["_N_int_medio"] = float(np.mean(n_int))
    return out


def fit(xs: list[float], ys: list[float]) -> tuple[float, float]:
    x, y = np.log(xs), np.log(ys)
    a = np.vstack([x, np.ones_like(x)]).T
    coef = np.linalg.lstsq(a, y, rcond=None)[0]
    r = y - a @ coef
    return float(coef[0]), math.sqrt(float(r @ r) / (len(x) - 2) * float(np.linalg.inv(a.T @ a)[0, 0]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["c4b", "small", "n2048", "final"], required=True)
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "final":
        c4b = json.loads((OUT / "c4b.json").read_text())["data"]
        d: dict[str, Any] = {}
        for p in ("small", "n2048"):
            d.update(json.loads((OUT / f"{p}.json").read_text())["data"])
        fits = {}
        for k in schedules(256):
            ei, si = fit([float(n) for n in NS], [d[str(n)][k]["D_int"] for n in NS])
            et, st = fit([float(n) for n in NS], [d[str(n)][k]["D_tot"] for n in NS])
            fits[k] = {"eta_int": ei, "se_int": si, "eta_tot": et, "se_tot": st, "D_int_2048": d["2048"][k]["D_int"],
                       "eta_int_compatible_con_0.5": abs(ei - 0.5) <= 2 * si}
        compat = [k for k, f in fits.items() if f["D_int_2048"] <= D_MAX and f["eta_int"] + 2 * f["se_int"] <= 0]
        incompat = all(f["eta_int"] - 2 * f["se_int"] > 0 for f in fits.values())
        verdict = ("INDETERMINADO (C4b no aprobado)" if not c4b["pass"]
                   else "COMPATIBLE" if compat else "INCOMPATIBLE" if incompat else "INDETERMINADO")
        summary = {"commit": head, "C4b": c4b["pass"], "fits": fits, "L4b": verdict}
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    data = part_c4b() if args.part == "c4b" else part_l4b((256, 512, 1024) if args.part == "small" else (2048,))
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s", data.get("pass", ""))


if __name__ == "__main__":
    main()
