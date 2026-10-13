"""E1″ (docs/OMEGA_E1DOBLEPRIMA_PREREGISTRO.md). Constantes congeladas.

Tramos:      --fam <familia> --n <N> --seeds <a> <b>
Evaluacion:  --final
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
from e1prima import criterion  # noqa: E402

from omega.causality.entropy import interval_balanced, knee_from_spectrum, pauli_jordan, sj_from_eig, ssee  # noqa: E402

EXPERIMENT_ID = 313
FAMS = {"geometrico": (20261054, 30), "kr": (20261055, 10), "percolacion": (20261056, 10), "capas_gruesas": (20261057, 10)}
NS = (800, 1200, 1600, 2000, 2400)
BETA_WIN = (-1.104, -0.904)
KNEE_WIN = (0.227, 0.906)
PREREG = Path("docs") / "OMEGA_E1DOBLEPRIMA_PREREGISTRO.md"
OUT = Path("results") / "omega_e1pp"


def one_run(fam: str, n: int, seed: int) -> dict[str, Any]:
    ent, _ = FAMS[fam]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ent, spawn_key=(EXPERIMENT_ID, n, seed))))
    c = FAMILIES[fam](n, rng)
    r = interval_balanced(c, n // 8)
    row: dict[str, Any] = {"fam": fam, "N": n, "seed": seed, "n_R": int(r.size)}
    mu, u = np.linalg.eigh(pauli_jordan(c))
    kn, beta = knee_from_spectrum(mu, n)
    row["rodilla"], row["beta"] = kn, beta
    if r.size == 0:
        return row
    idt, w = sj_from_eig(mu, u, math.sqrt(n) / (4 * math.pi))
    row["S_c1.0"] = ssee(idt, w, r, True, math.sqrt(r.size) / (4 * math.pi))
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fam")
    ap.add_argument("--n", type=int)
    ap.add_argument("--seeds", type=int, nargs=2)
    ap.add_argument("--final", action="store_true")
    args = ap.parse_args()
    if "C5 DECISIÓN: EJECUTAR" not in (Path("docs") / "OMEGA_CONSEJO_REV5_REPLICA.md").read_text():
        sys.exit("E1″ requiere la decisión del Consejo")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    if not args.final:
        t0 = time.time()
        rows = [one_run(args.fam, args.n, s) for s in range(args.seeds[0], args.seeds[1])]
        p = OUT / f"{args.fam}_{args.n}_{args.seeds[0]:02d}-{args.seeds[1]:02d}.json"
        p.write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "rows": rows}, default=str))
        print(p.name, "ok", round(time.time() - t0, 1), "s")
        return
    data: dict[str, list[dict[str, Any]]] = {}
    for fam, (_, seeds) in FAMS.items():
        data[fam] = [r for p in sorted(OUT.glob(f"{fam}_*.json")) for r in json.loads(p.read_text())["rows"]]
        if len(data[fam]) != seeds * len(NS):
            sys.exit(f"{fam}: {len(data[fam])} corridas de {seeds * len(NS)}")
    fits = {f: criterion(rows, 1.0) for f, rows in data.items()}
    g = fits["geometrico"]
    control = bool(g.get("r2_log", 0) > g.get("r2_lin", 1) and abs(g["a"] - 1 / 3) <= 2 * g["se"])
    match = {}
    for f in ("kr", "percolacion", "capas_gruesas"):
        ff = fits[f]
        if "a" not in ff:
            match[f] = False
            continue
        match[f] = bool(ff["r2_log"] > ff["r2_lin"] and abs(ff["a"] - g["a"]) <= 2 * math.sqrt(g["se"] ** 2 + ff["se"] ** 2))
    principal = "INDETERMINADO" if not control else ("NO DISCRIMINA" if any(match.values()) else "DISCRIMINA")
    spec = {}
    for f, rows in data.items():
        b = float(np.nanmedian([r["beta"] for r in rows]))
        ks = [r["rodilla"] for r in rows if r["rodilla"] is not None]
        k = float(np.median(ks)) if ks else None
        spec[f] = {"beta_mediana": b, "rodilla_mediana": k,
                   "compatible": bool(BETA_WIN[0] <= b <= BETA_WIN[1] and k is not None and KNEE_WIN[0] <= k <= KNEE_WIN[1])}
    secundario = ("INDETERMINADO" if not spec["geometrico"]["compatible"]
                  else "NO DISCRIMINA" if any(spec[f]["compatible"] for f in spec if f != "geometrico") else "DISCRIMINA")
    if principal == "INDETERMINADO":
        combinada = "INDETERMINADO"
    elif principal == "DISCRIMINA" and secundario == "DISCRIMINA":
        combinada = "DISCRIMINA (entropía y huella)"
    elif principal == "DISCRIMINA":
        combinada = "DISCRIMINA PARCIAL (solo entropía)"
    elif secundario == "DISCRIMINA":
        combinada = "DISCRIMINA PARCIAL (solo huella)"
    else:
        combinada = "NO DISCRIMINA"
    out = {"commit": head, "ajustes": fits, "control_valido": control, "coincide_con_geometria": match,
           "principal": principal, "huella": spec, "secundario": secundario, "E1_doble_prima": combinada}
    (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({k: out[k] for k in ("principal", "secundario", "E1_doble_prima")}, indent=2))


if __name__ == "__main__":
    main()
