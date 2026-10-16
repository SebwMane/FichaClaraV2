"""Omega-P3 L-P3-0b (docs/OMEGA_P3_L0.md seccion 6): curvatura de Ollivier-Ricci sobre un panel de grafos fijos, sin dinamica.

Panel de L-OD-0 (reutilizado de tools/omega_d_l0.py) sin la union de cliques, mas arbol aleatorio N=1000; clave PCG64
(20261011, familia, semilla). Por grafo: kappa perezoso (0.5, coste de saltos, componente gigante) en hasta 1500 aristas.
Criterios KD1-KD3 de la seccion 6. Los finales C0 (N = 729) son informativos.
Salida: results/p3_l0b/{graphs.jsonl, summary.json}. `--smoke`: 3 grafos pequenos -> results/p3_l0b_smoke/.
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
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import omega_d_l0 as D0  # noqa: E402
from c0_dynamics import _atomic_write  # noqa: E402
from p1_diagnostic import patchwork  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import strong_support  # noqa: E402
from omega.controls.random_geometric import rgg_torus  # noqa: E402
from omega.curvature.ollivier import ollivier_edge  # noqa: E402
from omega.diffusion.overlap import giant_dense  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402
from omega.geometry.distances import hop_distance_matrix  # noqa: E402
from omega.landscape.references import connected_caveman  # noqa: E402

MASTER = 20261011
N = 1000
IDLENESS = 0.5
MAX_EDGES = 1500
TAIL = 0.1
QUANTS = (10, 25, 75, 90)
TREE_FID = 15
OUT = ROOT / "results" / "p3_l0b"
OUT_SMOKE = ROOT / "results" / "p3_l0b_smoke"
_SPECS: list[dict[str, Any]] = []


def specs() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for s in D0.specs():
        if s["family"] == "cliques_K10":
            continue  # sustituida por caveman K8 (ya presente en el panel de L-OD-0)
        if s["key"] is not None:
            s = {**s, "key": (MASTER, s["key"][1], s["key"][2])}
        out.append(s)
    for sd in (0, 1, 2):
        out.append({"id": f"arbol_s{sd}", "family": "arbol", "cls": "NL", "gen": lambda g: R.random_tree(N, g), "key": (MASTER, TREE_FID, sd)})
    return out


def specs_smoke() -> list[dict[str, Any]]:
    n = 216
    return [
        {"id": "RGG3_k12_s0", "family": "RGG3_k12", "cls": "G", "gen": lambda g: rgg_torus(n, 3, 12.0, g), "key": (MASTER, 1, 0)},
        {"id": "retazos_2x3_s0", "family": "retazos_2x3", "cls": "A", "gen": lambda g: patchwork(n, 2, g), "key": (MASTER, 6, 0)},
        {"id": "caveman_K8_s0", "family": "caveman_K8", "cls": "R", "gen": lambda g: connected_caveman(200, 8), "key": (MASTER, 11, 0)},
    ]


def curvature_sample(adj: np.ndarray, key: tuple[int, ...] | None, max_edges: int = MAX_EDGES) -> dict[str, Any]:
    a = giant_dense(adj)
    n = int(a.shape[0])
    xs, ys = np.nonzero(np.triu(a, 1))
    m = int(xs.size)
    rng = R.rng_from_key((*key, 3) if key is not None else (MASTER, 99, 99, 3))
    if m > max_edges:
        pick = np.sort(rng.choice(m, size=max_edges, replace=False))
        xs, ys = xs[pick], ys[pick]
    w = a.astype(np.float64)
    d = hop_distance_matrix(a)
    kap = np.asarray([ollivier_edge(w, a, d, int(x), int(y), IDLENESS) for x, y in zip(xs, ys, strict=True)], dtype=np.float64)
    q = np.percentile(kap, QUANTS)
    return {"n_giant": n, "n_edges_total": m, "n_sampled": int(kap.size), "median_kappa": float(np.median(kap)),
            "mean_kappa": float(kap.mean()), "f_neg": float(np.mean(kap < -TAIL)), "f_pos": float(np.mean(kap > TAIL)),
            "quantiles": {str(p): float(v) for p, v in zip(QUANTS, q)}}


def run_idx(i: int) -> dict[str, Any]:
    spec = _SPECS[i]
    if "npz" in spec:
        z = np.load(spec["npz"])
        adj = strong_support(np.asarray(z[z.files[0]], dtype=np.float64)).astype(np.float64)
    else:
        adj = spec["gen"](R.rng_from_key(spec["key"]))
    t0 = time.time()
    row = {k: v for k, v in spec.items() if k not in ("gen", "npz")}
    row.update(curvature_sample(adj, spec["key"]))
    row["seconds"] = time.time() - t0
    return row


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by[r["family"]].append(r)
    med = {f: [r["median_kappa"] for r in rs] for f, rs in by.items()}
    res: dict[str, Any] = {"family_table": {f: {"n_graphs": len(rs), "n_giant": [r["n_giant"] for r in rs], "median_kappa": med[f],
                                                  "f_neg": [r["f_neg"] for r in rs], "f_pos": [r["f_pos"] for r in rs]}
                                            for f, rs in by.items()}}
    need1 = ("caveman_K8", "ER_k12", "RR_k12", "arbol", "T3_10", "cuadrado_32")
    if all(f in by for f in need1):
        cav = med["caveman_K8"]
        c_cav = all(v >= 0.2 for v in cav)
        c_nl = all(v <= -0.2 for f in ("ER_k12", "RR_k12", "arbol") for v in med[f])
        c_flat = all(abs(v) <= 0.02 for f in ("T3_10", "cuadrado_32") for v in med[f])
        res["KD1"] = {"pass": bool(c_cav and c_nl and c_flat), "caveman_ok": bool(c_cav), "ER_RR_arbol_ok": bool(c_nl),
                      "T3_cuadrado_ok": bool(c_flat), "medians": {f: med[f] for f in need1}}
    else:
        res["KD1"] = None
    rgg = by.get("RGG3_k12", [])
    if rgg and "retazos_2x3" in by and "retazos_3x3" in by:
        base = float(np.median([r["f_neg"] for r in rgg]))
        fam = {f: [r["f_neg"] for r in by[f]] for f in ("retazos_2x3", "retazos_3x3")}
        ok = all(len(v) == 3 and all(x >= base + 0.05 for x in v) for v in fam.values())
        res["KD2"] = {"pass": bool(ok), "f_neg_RGG3_median": base, "threshold": base + 0.05, "f_neg_retazos": fam}
    else:
        res["KD2"] = None
    if rgg and "RGG2_k12" in by:
        m3, m2 = float(np.median(med["RGG3_k12"])), float(np.median(med["RGG2_k12"]))
        res["KD3"] = {"informative": True, "median_RGG3": m3, "median_RGG2": m2, "RGG3_ge_0.1": bool(m3 >= 0.1),
                      "RGG2_ge_0.1": bool(m2 >= 0.1), "favors_crystals": bool(m3 >= 0.1 or m2 >= 0.1)}
    else:
        res["KD3"] = None
    k1, k2 = res["KD1"], res["KD2"]
    if k1 is None or k2 is None:
        dec = "INCOMPLETO (panel parcial)"
    elif not k1["pass"]:
        dec = "D muere"
    elif not k2["pass"]:
        dec = "D debilitada"
    else:
        dec = "premisa de D viable; proponer L-P3-1"
    res["decision"] = dec
    return res


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)
    out = OUT_SMOKE if args.smoke else OUT
    _SPECS[:] = specs_smoke() if args.smoke else specs()
    t0 = time.time()
    with mp.get_context("fork").Pool(args.procs) as pool:
        rows = pool.map(run_idx, range(len(_SPECS)), chunksize=1)
    out.mkdir(parents=True, exist_ok=True)
    _atomic_write(out / "graphs.jsonl", "".join(json.dumps(r, allow_nan=True) + "\n" for r in rows))
    res = evaluate([r for r in rows if r["cls"] != "INFO"])
    c0 = [r for r in rows if r["cls"] == "INFO"]
    res["C0_informative"] = {"n_graphs": len(c0), "ids": [r["id"] for r in c0], "median_kappa": [r["median_kappa"] for r in c0],
                             "f_neg": [r["f_neg"] for r in c0], "f_pos": [r["f_pos"] for r in c0]}
    summary = {"step": "p3_l0b", "smoke": args.smoke, "code_commit": head_commit(), "n_graphs": len(rows),
               "seconds": time.time() - t0, "result": res}
    _atomic_write(out / "summary.json", json.dumps(summary, indent=2, allow_nan=True) + "\n")
    for k in ("KD1", "KD2", "KD3"):
        print(k, None if res[k] is None else res[k].get("pass", res[k]))
    print("decision:", res["decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
