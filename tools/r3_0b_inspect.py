"""OMEGA R3-0b: inspeccion del residuo J3 (docs/OMEGA_R3_0.md §8, congelado). Reutiliza tools/r3_0b.py.

Salida: results/r3_0b_inspect/{runs.jsonl, summary.json}. Uso: python tools/r3_0b_inspect.py [--procs 3] [--fresh]
Reanudable por unidad (runs.jsonl).

Decisiones de implementacion
----------------------------
* I1: J3 semillas 3-11 con r3_0b.run_pair (par original de instancias independientes; claves (MASTER_R3, 8, semilla, 0/1)).
  Prediccion: std de delta >= 0.08 y media en [0.10, 0.25] sobre las semillas 3-11 (std muestral, ddof=1). Se informa
  ademas la estadistica con las 12 semillas (0-2 de results/r3_0b si existen) solo como dato.
* I2: semillas 0-2; UN solo poset por semilla, clave de generacion (MASTER_R3, 8, semilla, 0) para los tres tamanos
  (J3 es consistente por prefijos: el poset de n elementos es el inducido por los n primeros). Objetivos 1e4, 8e4, 6.4e5;
  tamano = menor n con |L| >= objetivo (|L| monotona en n por inclusion de subposets inducidos). Clave de MEDICION
  rc3.measure = (MASTER_R3, 8, semilla, k) con k = 0, 1, 2 el indice de tamano (solo distingue el RNG de la medida; el
  poset es el mismo). delta(N,8N) y delta(8N,64N) con rc3.scaling_delta sobre n_giant y R de half_mass_radius; exclusiones
  con pair_exclusions; W5 entre tamanos consecutivos (|E/N' - E/N|/(E/N) < 0.25).
  Prediccion I2: delta(N,8N) en [0.13,0.27] y delta(8N,64N) < delta(N,8N); se informa por semilla y "todas".
* I3: en lugar de usar solo ensayos de biseccion (truncados al objetivo, no dan |L| exacto) se calcula |L| EXACTO para
  cada n = 5, 6, ..., n(64N) con la misma clave de I2 (cada uno una enumeracion completa sin aristas); esos son los
  puntos (n, |L|) registrados. Ajustes MCO con intercepto: ln|L| = a + b sqrt(n) y ln|L| = a + b ln n; se informa SSE de
  ambos por semilla; "favorece sqrt" si SSE_sqrt < SSE_ln. Prediccion I3: sqrt favorecido (en las 3 semillas; se
  informa el recuento).
* Regla de clasificacion (§8): (a) delta(8N,64N) < delta(N,8N) en >= 2/3 semillas con descenso medio (sobre las 3) >= 0.02;
  (b) I3 favorece sqrt en >= 2/3 semillas. Ambas -> CRECIMIENTO INTERMEDIO; si no, NO CLASIFICADO.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import r3_0b as R  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402
from rc3 import measure, pair_exclusions, scaling_delta  # noqa: E402

from omega.experiments.v11.gate import head_commit  # noqa: E402

OUT = ROOT / "results" / "r3_0b_inspect"
FID_J3 = R.FAMILY_ID["J3"]
TARGETS3 = (10_000, 80_000, 640_000)
N_START = 5


def fit_sse(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """MCO y = a + b x; devuelve a, b, SSE."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.stack([np.ones_like(x), x], axis=1)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    return {"a": float(coef[0]), "b": float(coef[1]), "sse": float(res @ res)}


def growth_fits(points: list[tuple[int, int]]) -> dict[str, Any]:
    n = np.array([p[0] for p in points], float)
    y = np.log(np.array([p[1] for p in points], float))
    fs, fl = fit_sse(np.sqrt(n), y), fit_sse(np.log(n), y)
    return {"sqrt": fs, "log": fl, "favors_sqrt": bool(fs["sse"] < fl["sse"]), "n_points": len(points)}


def unit_i1(seed: int, cap: int) -> dict[str, Any]:
    spec = {"id": f"J3_s{seed}", "family": "J3", "gen": "J3", "param": 0, "seed": seed}
    row = R.run_pair(spec, cap)
    row["id"] = f"I1_s{seed}"
    row["unit"] = "I1"
    return row


def unit_i2(seed: int, cap: int) -> dict[str, Any]:
    t0 = time.perf_counter()
    key0 = (R.MASTER_R3, FID_J3, seed, 0)
    points: list[tuple[int, int]] = []
    n = N_START
    status = "OK"
    while True:
        en = R.enumerate_downsets(R.make_poset("J3", 0, n, key0), cap, edges=False)
        if not en["complete"]:
            status = "CAP_EXCEEDED"
            break
        points.append((n, en["count"]))
        if en["count"] >= TARGETS3[-1]:
            break
        n += 1
    t_pts = round(time.perf_counter() - t0, 1)
    row: dict[str, Any] = {"id": f"I2_s{seed}", "unit": "I2", "seed": seed, "status": status, "points": points,
                           "seconds_points": t_pts}
    if status != "OK":
        row["seconds"] = round(time.perf_counter() - t0, 1)
        return row
    sizes = []
    for k, tgt in enumerate(TARGETS3):
        n_k = next(nn for nn, L in points if L >= tgt)
        preds = R.make_poset("J3", 0, n_k, key0)
        en = R.enumerate_downsets(preds, cap)
        adj = R.cover_graph(en["count"], en["eu"], en["ev"])
        mkey = (R.MASTER_R3, FID_J3, seed, k)
        m = measure(adj, mkey)
        del adj
        m.update(key=list(mkey), n_elements=n_k, L=en["count"], target=tgt, E_over_N=m["mean_degree"] / 2.0)
        sizes.append(m)
    pairs = []
    for a, b in ((sizes[0], sizes[1]), (sizes[1], sizes[2])):
        d = scaling_delta(a["R"], b["R"], a["n_giant"] or 0, b["n_giant"] or 0)
        ex = pair_exclusions(d, (a["annulus_status"], b["annulus_status"]), (a["coherence_status"], b["coherence_status"]))
        pairs.append({"delta": d, "exclusions": ex,
                      "w5_valid": bool(abs(b["E_over_N"] - a["E_over_N"]) / a["E_over_N"] < 0.25)})
    row.update(sizes=sizes, pairs=pairs, delta1=pairs[0]["delta"], delta2=pairs[1]["delta"], fits=growth_fits(points))
    row["seconds"] = round(time.perf_counter() - t0, 1)
    return row


def _unit(u: tuple[str, int, int]) -> dict[str, Any]:
    kind, seed, cap = u
    return unit_i1(seed, cap) if kind == "I1" else unit_i2(seed, cap)


def _load(path: Path) -> dict[str, dict[str, Any]]:
    done: dict[str, dict[str, Any]] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                try:
                    r = json.loads(line)
                    done[r["id"]] = r
                except json.JSONDecodeError:
                    pass
    return done


def evaluate(rows: dict[str, dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    i1 = [rows[f"I1_s{s}"] for s in range(3, 12) if f"I1_s{s}" in rows]
    ds = [r["delta"] for r in i1 if r["delta"] is not None]
    if len(ds) >= 2:
        mean, std = float(np.mean(ds)), float(np.std(ds, ddof=1))
        out["I1"] = {"deltas": {r["id"]: r["delta"] for r in i1}, "n": len(ds), "mean": mean, "std": std,
                     "pred_std_ge_0.08": std >= 0.08, "pred_mean_in_[0.10,0.25]": 0.10 <= mean <= 0.25,
                     "match": bool(std >= 0.08 and 0.10 <= mean <= 0.25),
                     "excluded": [r["id"] for r in i1 if r["excluded"]], "failures": [r["id"] for r in i1 if r["status"] != "OK"]}
    i2 = [rows[f"I2_s{s}"] for s in range(3) if f"I2_s{s}" in rows and rows[f"I2_s{s}"]["status"] == "OK"]
    if i2:
        drops = [r["delta1"] - r["delta2"] for r in i2 if r["delta1"] is not None and r["delta2"] is not None]
        per = [{"id": r["id"], "n": [s["n_elements"] for s in r["sizes"]], "L": [s["L"] for s in r["sizes"]],
                "mean_degree": [s["mean_degree"] for s in r["sizes"]], "R": [s["R"] for s in r["sizes"]],
                "delta1": r["delta1"], "delta2": r["delta2"], "exclusions": [p["exclusions"] for p in r["pairs"]],
                "w5_valid": [p["w5_valid"] for p in r["pairs"]],
                "annulus": [s["annulus_status"] for s in r["sizes"]], "coherence": [s["coherence_status"] for s in r["sizes"]],
                "d1_in_[0.13,0.27]": r["delta1"] is not None and 0.13 <= r["delta1"] <= 0.27,
                "d2_lt_d1": r["delta1"] is not None and r["delta2"] is not None and r["delta2"] < r["delta1"]} for r in i2]
        n_dec = sum(p["d2_lt_d1"] for p in per)
        mean_drop = float(np.mean(drops)) if drops else None
        out["I2"] = {"per_seed": per, "n_decreasing": n_dec, "mean_drop": mean_drop,
                     "match_d1_range_all": all(p["d1_in_[0.13,0.27]"] for p in per) and len(per) == 3,
                     "match_decrease_all": n_dec == len(per) == 3}
        fav = [r["fits"]["favors_sqrt"] for r in i2]
        out["I3"] = {"per_seed": [{"id": r["id"], "n_points": r["fits"]["n_points"], "sse_sqrt": r["fits"]["sqrt"]["sse"],
                                   "sse_log": r["fits"]["log"]["sse"], "fit_sqrt": [r["fits"]["sqrt"]["a"], r["fits"]["sqrt"]["b"]],
                                   "fit_log": [r["fits"]["log"]["a"], r["fits"]["log"]["b"]],
                                   "favors_sqrt": r["fits"]["favors_sqrt"]} for r in i2],
                     "n_favor_sqrt": int(sum(fav)), "match_all": bool(all(fav) and len(fav) == 3)}
        a = bool(n_dec >= 2 and mean_drop is not None and mean_drop >= 0.02 and len(i2) == 3)
        b = bool(sum(fav) >= 2 and len(i2) == 3)
        out["classification"] = {"a": a, "b": b, "result": "CRECIMIENTO INTERMEDIO" if a and b else "NO CLASIFICADO"}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--cap", type=int, default=R.CAP)
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    pj = OUT / "runs.jsonl"
    if args.fresh and pj.exists():
        pj.unlink()
    done = _load(pj)
    units = [("I2", s, args.cap) for s in range(3)] + [("I1", s, args.cap) for s in range(3, 12)]
    todo = [u for u in units if f"{u[0]}_s{u[1]}" not in done]
    t0 = time.perf_counter()
    print(f"{len(units)} unidades ({len(done)} hechas), procs={args.procs}", flush=True)
    with mp.get_context("fork").Pool(args.procs) as pool, open(pj, "a") as fh:
        for row in pool.imap_unordered(_unit, todo, chunksize=1):
            fh.write(json.dumps(row, default=str) + "\n")
            fh.flush()
            done[row["id"]] = row
            if row["unit"] == "I1":
                print(f"  {row['id']} n={row['sizes'][0]['n_elements']}/{row['sizes'][1]['n_elements']} delta={row['delta']} "
                      f"{row['verdict']} W5={row['w5_valid']} {row['seconds']}s", flush=True)
            else:
                print(f"  {row['id']} {row['status']} d1={row.get('delta1')} d2={row.get('delta2')} {row['seconds']}s", flush=True)
    ev = evaluate(done)
    summary = {"step": "r3_0b_inspect", "code_commit": head_commit(), "seconds_this_run": round(time.perf_counter() - t0, 1),
               "evaluation": ev}
    _atomic_write(OUT / "summary.json", json.dumps(summary, indent=2, default=str))
    print(json.dumps(ev, indent=1, default=str)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
