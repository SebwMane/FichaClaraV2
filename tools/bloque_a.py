"""Bloque A (Rev. 3 §13 y §15): via analitica de la rama lorentziana. Constantes congeladas. No ajustar."""

from __future__ import annotations

import argparse
import cmath
import itertools
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import between_counts, bd_action_2d_smeared, smearing_weight  # noqa: E402
from omega.causality.orders import random_2order, random_bipartite, transitive_closure  # noqa: E402
from omega.causality.phase import phase_average  # noqa: E402

EXPERIMENT_ID = 305
ENT_A2, ENT_A3, ENT_GAUSS = 20261029, 20261030, 20261031
OUT = Path("results") / "bloque_a"


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def smeared_from_between(vals: np.ndarray[Any, Any], n: int, eps: float) -> float:
    return float(eps * (n - 2 * eps * float(np.sum(smearing_weight(vals.astype(np.float64), eps)))))


def a2a(m: int) -> dict[str, Any]:
    rows = []
    for i, n in enumerate((8, 12)):
        orders = [random_bipartite(n, 0.5, rng_for(ENT_A2, i, s)) for s in range(m)]
        for eps in (0.5, 0.25):
            acts = [bd_action_2d_smeared(c, eps) for c in orders]
            for x in (0.1, 0.3, 0.6):
                lam = x / eps**2
                phi, se = phase_average(acts, lam)
                exact = abs(math.cos(x)) ** ((n // 2) * (n - n // 2))
                rows.append({"n": n, "eps": eps, "lam_eps2": x, "phi": phi, "se": se, "exact": exact, "ok": abs(phi - exact) <= 3 * se})
    return {"rows": rows, "pass": all(r["ok"] for r in rows)}


def a2b() -> dict[str, Any]:
    n, l1, l2 = 8, [0, 1], [2, 3, 4, 5]
    l3 = [6, 7]
    pairs = [(a, b) for a in l1 for b in l2] + [(a, b) for a in l2 for b in l3]
    epss = (0.5, 0.25, 0.1)
    xs = (0.3, math.pi / 3 + 0.1, math.pi / 2, 2 * math.pi / 3 - 0.1)
    acts: dict[float, list[float]] = {e: [] for e in epss}
    links: list[int] = []
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        adj = np.zeros((n, n), dtype=bool)
        for (a, b), on in zip(pairs, bits):
            adj[a, b] = bool(on)
        c = transitive_closure(adj)
        bvals = between_counts(c)
        vals = bvals[bvals >= 0]
        for e in epss:
            acts[e].append(smeared_from_between(vals, n, e))
        links.append(sum(bits))
    lk = np.array(links, dtype=np.float64)
    rows = []
    for e in epss:
        s = np.array(acts[e])
        delta = s - (e * n - 2 * e * e * lk)
        for x in xs:
            lam = x / e**2
            phi_enum = complex(np.mean(np.exp(1j * lam * s)))
            phi_lead = cmath.exp(1j * lam * e * n) * ((1 + cmath.exp(-2j * x)) / 2) ** len(pairs)
            bound = float(np.mean(np.abs(lam * delta)))
            diff = abs(phi_enum - phi_lead)
            rows.append({"eps": e, "lam_eps2": x, "abs_enum": abs(phi_enum), "abs_lead": abs(phi_lead),
                         "diff": diff, "bound": bound, "ok": diff <= bound + 1e-12})
    return {"rows": rows, "pass": all(r["ok"] for r in rows), "configs": 2 ** len(pairs)}


def a3_exponent(m: int, ns: tuple[int, ...]) -> dict[str, Any]:
    var: dict[float, dict[int, float]] = {0.25: {}, 0.1: {}}
    for i, n in enumerate(ns):
        sa: dict[float, list[float]] = {0.25: [], 0.1: []}
        for s in range(m):
            b = between_counts(random_2order(n, rng_for(ENT_A3, i, s)))
            vals = b[b >= 0]
            for e in sa:
                sa[e].append(smeared_from_between(vals, n, e) / e**2)
        for e in sa:
            var[e][n] = float(np.var(sa[e], ddof=1))
    fits = {}
    for e, d in var.items():
        x = np.log(list(d.keys()))
        y = np.log(list(d.values()))
        a = np.vstack([x, np.ones_like(x)]).T
        coef = np.linalg.lstsq(a, y, rcond=None)[0]
        resid = y - a @ coef
        s2 = float(resid @ resid) / (len(x) - 2)
        se = math.sqrt(s2 * float(np.linalg.inv(a.T @ a)[0, 0]))
        fits[str(e)] = {"gamma": float(coef[0]), "se": se, "c": float(coef[1]), "var": {str(k): v for k, v in d.items()},
                        "ok": float(coef[0]) + 2 * se < 2}
    return fits


def a3_gauss(m: int) -> dict[str, Any]:
    rows = []
    thr = 3 / math.sqrt(m)
    for i, n in enumerate((16, 32, 64)):
        sa: dict[float, list[float]] = {0.25: [], 0.1: []}
        for s in range(m):
            b = between_counts(random_2order(n, rng_for(ENT_GAUSS, i, s)))
            vals = b[b >= 0]
            for e in sa:
                sa[e].append(smeared_from_between(vals, n, e))
        for e, acts in sa.items():
            v = float(np.var(np.array(acts) / e**2, ddof=1))
            se_v = v * math.sqrt(2 / (m - 1))
            for x in (0.02, 0.05, 0.1):
                lam = x / e**2
                phi, se = phase_average(acts, lam)
                pred = math.exp(-x * x * v / 2)
                tol = 3 * se + 3 * abs(-x * x / 2 * pred) * se_v
                usable = phi > thr
                rows.append({"n": n, "eps": e, "lam_eps2": x, "phi": phi, "pred": pred, "tol": tol, "usable": usable,
                             "ok": (not usable) or abs(phi - pred) <= tol})
    return {"rows": rows, "pass": all(r["ok"] for r in rows), "n_usable": sum(r["usable"] for r in rows)}


def crossover(fits: dict[str, Any]) -> dict[str, Any]:
    x = 1.2
    out = {}
    for e, f in fits.items():
        nstar = None
        for n in range(16, 200001, 16):
            geo = math.lgamma(n + 1) - x * x * math.exp(f["c"]) * n ** f["gamma"] / 2
            lay = (n * n / 4) * math.log(abs(2 * math.cos(x)))
            if geo > lay:
                nstar = n
                break
        out[e] = nstar
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--part", choices=["a2", "a3", "final"], default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "smoke":
        t0 = time.time()
        r = {"a2a": a2a(500)["pass"], "a3": {k: (round(v["gamma"], 2), round(v["se"], 2)) for k, v in a3_exponent(20, (64, 128, 256)).items()},
             "gauss_usable": a3_gauss(500)["n_usable"], "s": round(time.time() - t0, 1)}
        print(r)
        return
    if dirty:
        sys.exit("modo full exige arbol git limpio (fuera de results/)")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.part == "final":
        a2 = json.loads((OUT / "a2.json").read_text())
        a3 = json.loads((OUT / "a3.json").read_text())
        exp_ok = all(f["ok"] for f in a3["data"]["exponent"].values())
        summary = {
            "phase": "final", "commits": sorted({a2["commit"], a3["commit"]}),
            "A-2a": a2["data"]["a2a"]["pass"], "A-2b": a2["data"]["a2b"]["pass"],
            "A-3_exponente": {k: (v["gamma"], v["se"]) for k, v in a3["data"]["exponent"].items()},
            "A-3_gauss": a3["data"]["gauss"]["pass"], "A-3_gauss_puntos_utilizables": a3["data"]["gauss"]["n_usable"],
            "A-3": "PASA" if exp_ok and a3["data"]["gauss"]["pass"] else "NO PASA",
            "N_estrella_nivel4": crossover(a3["data"]["exponent"]),
        }
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    t0 = time.time()
    if args.part == "a2":
        data: dict[str, Any] = {"a2a": a2a(20000), "a2b": a2b()}
    else:
        data = {"exponent": a3_exponent(200, (64, 128, 256, 512, 1024)), "gauss": a3_gauss(20000)}
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
