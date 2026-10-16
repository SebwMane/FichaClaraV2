"""L-2c (Rev. 3 §10): compuerta de paisaje de Omega-C, ramas euclidea (E) y lorentziana (L), criterios separados.

Constantes congeladas en docs/OMEGA_2_PROPUESTA_CAUSAL.md §10. No ajustar.
Modos: smoke (rapido, sin veredicto) y full (arbol limpio, veredictos calculados).
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

from omega.causality.estimators import bd_action_2d, interval_abundances, mm_dimension  # noqa: E402
from omega.causality.orders import layered_2order, random_2order, random_bipartite  # noqa: E402
from omega.causality.phase import bipartite_phase_exact, phase_average  # noqa: E402
from omega.causality.profile import K_MAX, d_ref, profile_distance, reference_abundances  # noqa: E402

EXPERIMENT_ID = 302
ENT_E, ENT_LI, ENT_LIII = 20261018, 20261020, 20261019
_C04B = Path(__file__).resolve().parents[1] / "results" / "c0_4b_validation" / "summary.json"
_D2 = float(json.loads(_C04B.read_text())["delta"]["2"])  # AUDITORIA A-9: se lee del resultado oficial, no se copia
DELTA_C04B = {1: _D2, 2: _D2}
LAYERS = (2, 4, 16, 128)
LI_N, LI_LAM, LI_M = (8, 12), (0.2, 0.5, 1.0, math.pi / 2), 20000
LIII_N, LIII_LAM, LIII_M = (16, 32, 64, 128, 256), (1.2, math.pi / 2, 1.9), 20000


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def c04b_flag_2d(c: np.ndarray[Any, Any], ref2: np.ndarray[Any, Any]) -> bool:
    """Bandera C0-4b para ordenes con d_ref <= 2 (los unicos que aparecen en L-2c-E)."""
    dr = d_ref(mm_dimension(c))
    if dr is None or dr > 2:
        return False
    ref = ref2 if dr == 2 else reference_abundances(1, c.shape[0], K_MAX)
    return bool(profile_distance(interval_abundances(c, K_MAX), ref) <= DELTA_C04B[dr])


def branch_e(n: int, seeds: int) -> dict[str, Any]:
    ref2 = reference_abundances(2, n, K_MAX)
    geo: list[int] = []
    geo_flags: list[bool] = []
    for s in range(seeds):
        c = random_2order(n, rng_for(ENT_E, 1, s))
        f = c04b_flag_2d(c, ref2)
        geo_flags.append(f)
        if f:
            geo.append(bd_action_2d(c))
    nongeo: dict[str, Any] = {}
    for k in LAYERS:
        c = layered_2order(n, k, rng_for(ENT_E, 10 + k, 0))
        nongeo[f"k{k}"] = {"S": bd_action_2d(c), "flag_c04b": c04b_flag_2d(c, ref2)}
    usable = [v["S"] for v in nongeo.values() if not v["flag_c04b"]]
    delta = float(min(usable) - np.mean(geo)) if usable and geo else float("nan")
    return {"geo_S": geo, "geo_flag_fraction": float(np.mean(geo_flags)), "nongeo": nongeo, "delta_S_E": delta}


def branch_l_i() -> dict[str, Any]:
    rows = []
    for i, n in enumerate(LI_N):
        acts = [bd_action_2d(random_bipartite(n, 0.5, rng_for(ENT_LI, i, s))) for s in range(LI_M)]
        for lam in LI_LAM:
            mag, se = phase_average(acts, lam)
            exact = bipartite_phase_exact(n, lam)
            rows.append({"n": n, "lam": lam, "phi": mag, "se": se, "exact": exact, "ok": abs(mag - exact) <= 3 * se})
    return {"rows": rows, "pass": all(r["ok"] for r in rows)}


def branch_l_iii(ns: tuple[int, ...], m: int) -> dict[str, Any]:
    acts: dict[int, list[int]] = {}
    for i, n in enumerate(ns):
        acts[n] = [bd_action_2d(random_2order(n, rng_for(ENT_LIII, i, s))) for s in range(m)]
    fits: dict[str, Any] = {}
    thr = 3 / math.sqrt(m)
    for lam in LIII_LAM:
        pts = []
        for n in ns:
            mag, se = phase_average(acts[n], lam)
            pts.append({"n": n, "phi": mag, "se": se, "usable": mag > thr})
        # AUDITORIA A-10 (enmienda L2c-A1): phi >= 1 no admite -ln(phi) > 0 y se trata como no utilizable
        for p in pts:
            if p["phi"] >= 1 - 1e-12:
                p["usable"] = False
        use = [p for p in pts if p["usable"]]
        fit: dict[str, Any] = {"points": pts, "n_usable": len(use)}
        if len(use) >= 3:
            x = np.log([p["n"] for p in use])
            y = np.log([-math.log(p["phi"]) for p in use])
            a = np.vstack([x, np.ones_like(x)]).T
            coef, res, _, _ = np.linalg.lstsq(a, y, rcond=None)
            dof = len(use) - 2
            resid = y - a @ coef
            s2 = float(resid @ resid) / dof if dof > 0 else float("nan")
            cov = s2 * np.linalg.inv(a.T @ a)
            fit.update({"b": float(coef[0]), "b_se": float(math.sqrt(cov[0, 0])) if dof > 0 else float("nan")})
        fits[f"{lam:.4f}"] = fit
    return {"fits": fits, "action_sd": {str(n): float(np.std(acts[n], ddof=1)) for n in ns}}


def verdict_l_iii(res: dict[str, Any]) -> str:
    usable = [f for f in res["fits"].values() if f["n_usable"] >= 3]
    if not usable:
        return "INDETERMINADO"
    ok = [f for f in usable if math.isfinite(f.get("b_se", float("nan"))) and f["b"] + 2 * f["b_se"] < 2]
    return "PASA" if ok else "NO PASA"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "full" and dirty:
        sys.exit("modo full exige arbol git limpio")
    out = args.out or Path("results") / ("l2c" if args.mode == "full" else "l2c_smoke")
    out.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {"phase": "preregistered", "mode": args.mode, "commit": head, "dirty": dirty}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    t0 = time.time()
    if args.mode == "full":
        e = branch_e(1024, 20)
        li = branch_l_i()
        liii = branch_l_iii(LIII_N, LIII_M)
    else:
        e = branch_e(128, 4)
        li = {"rows": [], "pass": None}
        liii = branch_l_iii((16, 32, 64), 2000)
    with (out / "runs.jsonl").open("w") as fh:
        fh.write(json.dumps({"E": e, "L_i": li, "L_iii": liii}) + "\n")
    summary.update({"phase": "final", "seconds": round(time.time() - t0, 1)})
    if args.mode == "full":
        d = e["delta_S_E"]
        summary["verdicts"] = {
            # AUDITORIA A-11: un delta nan (sin geometricos marcados o sin referencias validas) era leido como "NO preferida"
            "L-2c-E": "INDETERMINADO" if not math.isfinite(d) else (
                "geometria preferida a beta->inf" if d > 0 else "geometria NO preferida a beta->inf"
            ),
            "delta_S_E": d,
            "L-i": li["pass"],
            "L-ii_ventana_KR": "(pi/3, 2pi/3) analitica",
            "L-iii": verdict_l_iii(liii),
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
