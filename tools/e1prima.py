"""E1′ (docs/OMEGA_E1PRIMA_PREREGISTRO.md). Constantes congeladas.

Uso por tramos:  --stage vprima|e1 --fam <familia> --n <N> --seeds <a> <b>
Evaluacion:      --final vprima|e1
La etapa e1 se niega a correr si el preregistro no contiene la linea "R4-1 FIRMADA".
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

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e0_e1 import FAMILIES  # noqa: E402

from omega.causality.entropy import (  # noqa: E402
    causal_complement,
    interval_balanced,
    knee_from_spectrum,
    pauli_jordan,
    sj_from_eig,
    ssee,
)

EXPERIMENT_ID = 312
ENT = {"vprima": 20261052, "e1": 20261053}
NS = (800, 1200, 1600, 2000, 2400)
SEEDS = 10
CS = (0.5, 1.0, 2.0)
SHAPE = (0.085, 0.125)
PREREG = Path("docs") / "OMEGA_E1PRIMA_PREREGISTRO.md"
OUT = Path("results") / "omega_e1p"


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def one_run(fam: str, n: int, seed: int, stage: str) -> dict[str, Any]:
    k = list(FAMILIES).index(fam)
    c = FAMILIES[fam](n, rng_for(ENT[stage], k, n, seed))
    r = interval_balanced(c, n // 8)
    row: dict[str, Any] = {"fam": fam, "N": n, "seed": seed, "n_R": int(r.size)}
    if r.size == 0:
        return row
    past = c.sum(axis=0)
    fut = c.sum(axis=1)
    row["frac_pasado"] = float(past[int(r[c[r][:, r].sum(axis=0).argmin()])] / n)
    row["frac_futuro"] = float(fut[int(r[c[r][:, r].sum(axis=1).argmin()])] / n)
    mu, u = np.linalg.eigh(pauli_jordan(c))
    for cc in CS:
        idt, w = sj_from_eig(mu, u, cc * math.sqrt(n) / (4 * math.pi))
        row[f"S_c{cc}"] = ssee(idt, w, r, True, cc * math.sqrt(r.size) / (4 * math.pi))
        if cc == 1.0:
            rp = causal_complement(c, r)
            row["n_Rperp"] = int(rp.size)
            row["S_perp_c1"] = ssee(idt, w, rp, True, math.sqrt(rp.size) / (4 * math.pi)) if rp.size else None
    kn, beta = knee_from_spectrum(mu, n)
    row["rodilla"], row["beta"] = kn, beta
    return row


def criterion(rows: list[dict[str, Any]], cc: float) -> dict[str, Any]:
    ok = [r for r in rows if r["n_R"] > 0]
    if len(ok) < 3:
        return {"pasa": False, "nota": "insuficientes corridas con region"}
    x = np.log(np.sqrt(np.array([r["n_R"] for r in ok], dtype=float)) / (4 * math.pi))
    xl = np.array([r["n_R"] for r in ok], dtype=float)
    y = np.array([r[f"S_c{cc}"] for r in ok], dtype=float)
    tot = float(np.sum((y - y.mean()) ** 2)) or 1e-300

    def fit(xx: np.ndarray[Any, Any]) -> tuple[float, float, float]:
        a = np.vstack([xx, np.ones_like(xx)]).T
        coef = np.linalg.lstsq(a, y, rcond=None)[0]
        res = y - a @ coef
        se = math.sqrt(float(res @ res) / (len(y) - 2) * float(np.linalg.inv(a.T @ a)[0, 0]))
        return float(coef[0]), se, 1 - float(res @ res) / tot

    a_log, se_log, r2_log = fit(x)
    _, _, r2_lin = fit(xl)
    pasa = r2_log > r2_lin and abs(a_log - 1 / 3) <= 2 * se_log + 0.05
    return {"pasa": bool(pasa), "a": a_log, "se": se_log, "r2_log": r2_log, "r2_lin": r2_lin}


def load(stage: str, fam: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for p in sorted(OUT.glob(f"{stage}_{fam}_*.json")):
        rows += json.loads(p.read_text())["rows"]
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["vprima", "e1"])
    ap.add_argument("--fam", default="geometrico")
    ap.add_argument("--n", type=int)
    ap.add_argument("--seeds", type=int, nargs=2)
    ap.add_argument("--final", choices=["vprima", "e1"])
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.final == "vprima":
        rows = load("vprima", "geometrico")
        if len(rows) != len(NS) * SEEDS:
            sys.exit(f"faltan corridas: {len(rows)} de {len(NS) * SEEDS}")
        shape_ok = float(np.mean([SHAPE[0] <= r.get("frac_pasado", -1) <= SHAPE[1] and SHAPE[0] <= r.get("frac_futuro", -1) <= SHAPE[1] for r in rows]))
        crit = criterion(rows, 1.0)
        out = {"commit": head, "forma_fraccion": shape_ok, "criterio_c1": crit, "V_prima": "PASA" if shape_ok >= 0.9 and crit["pasa"] else "FALLA"}
        (OUT / "vprima_summary.json").write_text(json.dumps(out, indent=2))
        print(json.dumps(out, indent=2))
        return
    if args.final == "e1" or args.stage == "e1":
        if "R4-1 FIRMADA" not in PREREG.read_text():
            sys.exit("E1′ requiere la firma R4-1 del usuario en el preregistro")
    if args.final == "e1":
        res: dict[str, Any] = {}
        for fam in FAMILIES:
            rows = load("e1", fam)
            if len(rows) != len(NS) * SEEDS:
                sys.exit(f"{fam}: faltan corridas ({len(rows)})")
            if all(r["n_R"] == 0 for r in rows if r["N"] == NS[-1]):
                res[fam] = {"ley_area": False, "nota": "no aplicable: sin regiones extensas"}
                continue
            per_c = {str(cc): criterion(rows, cc) for cc in CS}
            comp = [abs(r["S_c1.0"] - r["S_perp_c1"]) / r["S_c1.0"] for r in rows if r.get("S_perp_c1") and r.get("S_c1.0")]
            res[fam] = {"ley_area": per_c["1.0"]["pasa"], "barrido": per_c,
                        "complementariedad_media": float(np.mean(comp)) if comp else None,
                        "rodilla_mediana": float(np.median([r["rodilla"] for r in rows if r.get("rodilla")])) if any(r.get("rodilla") for r in rows) else None,
                        "beta_mediana": float(np.nanmedian([r["beta"] for r in rows if "beta" in r]))}
        geo = res["geometrico"]["ley_area"]
        others = [f for f in res if f != "geometrico" and res[f]["ley_area"]]
        verdict = "INDETERMINADO" if not geo else ("NO DISCRIMINA" if others else "DISCRIMINA")
        out = {"commit": head, "familias": res, "E1_prima": verdict}
        (OUT / "e1_summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps(out, indent=2, default=str))
        return
    t0 = time.time()
    rows = [one_run(args.fam, args.n, s, args.stage) for s in range(args.seeds[0], args.seeds[1])]
    path = OUT / f"{args.stage}_{args.fam}_{args.n}_{args.seeds[0]:02d}-{args.seeds[1]:02d}.json"
    path.write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "rows": rows}, default=str))
    print(path.name, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
